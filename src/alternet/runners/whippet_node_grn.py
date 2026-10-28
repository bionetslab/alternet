import argparse
import os.path as op

import pandas as pd
from distributed import Client, LocalCluster

from alternet.data_preprocessing import standardize_dataframe
from alternet.genome_annotation import read_annotation
from alternet.inference import compute_grn
from alternet.informative_exons import (informative_nodes, is_expressed, node_psi_matrix, read_exon_nodes,
                                        read_transcript_tpm, transcript_nodes, variable_nodes)


def psi_sample_name(psi_file):
    """Sample name of a Whippet output file <sample>.psi(.gz)."""
    return op.basename(psi_file).removesuffix('.gz').removesuffix('.psi')


def whippet_node_grn(psi_files, tpm_file, gtf_file, exon_nodes_file, regulator_file,
                     regulator_types=('TF', 'SF', 'TF_SF'), min_tpm=1.0, min_fraction=0.25, min_reads=5,
                     psi_low=0.05, psi_high=0.95, max_missing=0.5, num_permutations=1000, num_target_clusters=1,
                     seed=11, client=None):
    """
    GRN with the informative Whippet nodes of regulator genes (Psi) as regulators and the expressed
    transcripts (TPM) as targets, with signifikante P-values.

    Parameters:
        psi_files (list): Whippet <sample>.psi(.gz) files, one per sample. The sample names must match
            the columns of tpm_file.
        tpm_file (str): Merged Salmon transcript TPM table (see read_transcript_tpm).
        gtf_file (str): GTF used for the Whippet index.
        exon_nodes_file (str): Whippet exon-to-node mapping (graph.exons.tab(.gz)).
        regulator_file (str): Regulator list (data/regulator_db.tsv).
        regulator_types (list): Regulator_type values whose genes provide regulator nodes.
        min_tpm (float): Minimal TPM of an expressed transcript.
        min_fraction (float): Minimal fraction of samples in which a transcript is expressed / a node is variable.
        min_reads (int): Minimal Total_Reads to use the Psi of a node.
        psi_low (float): Psi below this counts as excluded.
        psi_high (float): Psi above this counts as included.
        max_missing (float): Maximal fraction of samples without a Psi value for a regulator node.
            Remaining missing values are filled with the node's mean Psi.
        num_permutations (int): Number of permutations for the empirical P-values.
        num_target_clusters (int): Number of target clusters for the permutation background.
        seed (int): Random seed for the regressors.
        client (dask.distributed.Client): Dask client. A new local cluster is started if None.

    Returns:
        tuple: (network, nodes). network has the columns of compute_grn, sources are named
            '<gene>:<node>'. nodes are the informative nodes (see informative_nodes) with their
            name in 'source'.
    """
    psi_files = {psi_sample_name(f): f for f in psi_files}
    tpm = read_transcript_tpm(tpm_file)
    samples = [s for s in tpm.columns if s in psi_files]
    if not samples:
        raise ValueError(f"No sample of the Psi files is a column of {tpm_file}")
    print(f"{len(samples)} samples with Psi and TPM ({len(psi_files)} Psi files, {tpm.shape[1]} TPM columns)")
    psi_files = {s: psi_files[s] for s in samples}
    tpm = tpm[samples]

    regulators = pd.read_csv(regulator_file, sep='\t')
    regulator_genes = set(regulators.loc[regulators['Regulator_type'].isin(regulator_types), 'Gene stable ID'])

    expressed = is_expressed(tpm, min_tpm, min_fraction).to_frame('study')
    variable = variable_nodes({'study': list(psi_files.values())}, min_fraction, min_reads, psi_low, psi_high)
    transcripts = transcript_nodes(read_annotation(gtf_file, features={'exon'}), read_exon_nodes(exon_nodes_file))
    nodes = informative_nodes(variable, expressed, transcripts, regulator_genes)
    nodes['source'] = nodes['Gene'].str.split('.').str[0] + ':' + nodes['Node'].astype(str)

    psi = node_psi_matrix(psi_files, nodes, min_reads)
    psi = psi.loc[:, psi.isna().mean() <= max_missing]
    psi = psi.fillna(psi.mean())
    psi.columns = [f"{gene.split('.')[0]}:{node}" for gene, node in psi.columns]
    print(f"{nodes.shape[0]} informative nodes, {psi.shape[1]} with Psi in enough samples")

    targets = tpm.loc[expressed['study']].T
    print(f"{targets.shape[1]} expressed target transcripts")
    data = standardize_dataframe(pd.concat([psi, targets], axis=1))

    network = compute_grn(data, list(targets.columns), list(psi.columns), client=client, seed=seed,
                          num_target_cluster=num_target_clusters, num_permutations=num_permutations)
    return network, nodes


def main():
    parser = argparse.ArgumentParser(
        description="GRN with the informative Whippet nodes of regulator genes (Psi) as regulators and the "
                    "expressed transcripts (TPM) as targets, with signifikante P-values."
    )
    parser.add_argument("--psi_files", type=str, nargs='+', required=True,
                        help="Whippet <sample>.psi(.gz) files, one per sample. The sample names must match "
                             "the columns of --tpm_file.")
    parser.add_argument("--tpm_file", type=str, required=True,
                        help="Merged Salmon transcript TPM table (comma-separated, transcript IDs in 'Name').")
    parser.add_argument("--gtf", type=str, required=True, help="GTF used for the Whippet index.")
    parser.add_argument("--exon_nodes", type=str, required=True,
                        help="Whippet exon-to-node mapping (graph.exons.tab(.gz) of the Whippet index).")
    parser.add_argument("--regulator_file", type=str, required=True, help="Regulator list (data/regulator_db.tsv).")
    parser.add_argument("--regulator_types", type=str, nargs='+', default=['TF', 'SF', 'TF_SF'],
                        help="Regulator_type values whose genes provide regulator nodes.")
    parser.add_argument("--output", type=str, required=True, help="Output network TSV.")
    parser.add_argument("--nodes_output", type=str, required=False,
                        help="Optional TSV of the informative nodes used as regulators.")
    parser.add_argument("--min_tpm", type=float, default=1.0, help="Minimal TPM of an expressed transcript.")
    parser.add_argument("--min_fraction", type=float, default=0.25,
                        help="Minimal fraction of samples in which a transcript is expressed / a node is variable.")
    parser.add_argument("--min_reads", type=int, default=5, help="Minimal Total_Reads to use the Psi of a node.")
    parser.add_argument("--psi_low", type=float, default=0.05, help="Psi below this counts as excluded.")
    parser.add_argument("--psi_high", type=float, default=0.95, help="Psi above this counts as included.")
    parser.add_argument("--max_missing", type=float, default=0.5,
                        help="Maximal fraction of samples without a Psi value for a regulator node. "
                             "Remaining missing values are filled with the node's mean Psi.")
    parser.add_argument("--num_permutations", type=int, default=1000)
    parser.add_argument("--num_target_clusters", type=int, default=1)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    client = Client(LocalCluster())
    network, nodes = whippet_node_grn(args.psi_files, args.tpm_file, args.gtf, args.exon_nodes, args.regulator_file,
                                      args.regulator_types, args.min_tpm, args.min_fraction, args.min_reads,
                                      args.psi_low, args.psi_high, args.max_missing, args.num_permutations,
                                      args.num_target_clusters, args.seed, client)
    client.close()

    network.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(network)} edges to {args.output}")
    if args.nodes_output is not None:
        nodes.to_csv(args.nodes_output, sep='\t', index=False)


if __name__ == "__main__":
    main()
