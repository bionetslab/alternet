import argparse

import pandas as pd
from distributed import Client, LocalCluster

from alternet.benchmark import Benchmark
from alternet.data_preprocessing import standardize_dataframe
from alternet.inference import compute_grn
from alternet.informative_exons import is_expressed, read_transcript_tpm


def gene_count_grn(tpm_file, biomart_file, regulator_file, regulator_types=('TF', 'SF', 'TF_SF'), min_tpm=1.0,
                   min_fraction=0.25, num_permutations=1000, num_target_clusters=1, seed=11, client=None,
                   benchmark=None, tpm_sep=','):
    """
    Gene GRN with regulator genes as regulators, with signifikante P-values. The expression of a
    gene is the sum of the TPM of its expressed transcripts.

    Parameters:
        tpm_file (str): Merged Salmon transcript TPM table (see read_transcript_tpm).
        biomart_file (str): Biomart table with 'Gene stable ID' and 'Transcript stable ID' (data/biomart.txt).
            Transcripts that are not in it are dropped.
        regulator_file (str): Regulator list (data/regulator_db.tsv).
        regulator_types (list): Regulator_type values whose genes are regulators.
        min_tpm (float): Minimal TPM of an expressed transcript.
        min_fraction (float): Minimal fraction of samples in which a transcript is expressed.
        num_permutations (int): Number of permutations for the empirical P-values.
        num_target_clusters (int): Number of target clusters for the permutation background.
        seed (int): Random seed for the regressors.
        client (dask.distributed.Client): Dask client. A new local cluster is started if None.
        benchmark (Benchmark): Records run time and memory of each step. A new one is used if None.
        tpm_sep (str): Column separator of tpm_file.

    Returns:
        pd.DataFrame: Network with the columns of compute_grn.
    """
    if benchmark is None:
        benchmark = Benchmark()

    with benchmark.step('read_inputs') as info:
        tpm = read_transcript_tpm(tpm_file, tpm_sep)
        tpm = tpm.loc[is_expressed(tpm, min_tpm, min_fraction)]
        biomart = pd.read_csv(biomart_file, sep='\t')
        tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
        genes = tpm.index.map(tx2gene)
        print(f"{tpm.shape[1]} samples, {tpm.shape[0]} expressed transcripts "
              f"({genes.isna().sum()} without a gene in {biomart_file})")
        gene_tpm = tpm.loc[genes.notna()].groupby(genes[genes.notna()]).sum().T
        print(f"{gene_tpm.shape[1]} genes with expressed transcripts")
        data = standardize_dataframe(gene_tpm)

        regulators = pd.read_csv(regulator_file, sep='\t')
        regulator_genes = list(set(regulators.loc[regulators['Regulator_type'].isin(regulator_types),
                                                  'Gene stable ID']) & set(data.columns))
        print(f"{len(regulator_genes)} expressed regulator genes")
        info.update(n_samples=data.shape[0], n_transcripts=tpm.shape[0], n_genes=data.shape[1])

    with benchmark.step('inference', n_regulators=len(regulator_genes), n_targets=data.shape[1],
                        num_permutations=num_permutations) as info:
        network = compute_grn(data, 'all', regulator_genes, client=client, seed=seed,
                              num_target_cluster=num_target_clusters, num_permutations=num_permutations)
        info.update(n_edges=len(network))
    return network


def main():
    parser = argparse.ArgumentParser(
        description="Gene GRN with regulator genes as regulators, with signifikante P-values. The expression "
                    "of a gene is the sum of the TPM of its expressed transcripts."
    )
    parser.add_argument("--tpm_file", type=str, required=True,
                        help="Merged Salmon transcript TPM table (transcript IDs in 'Name').")
    parser.add_argument("--tpm_sep", type=str, default=',',
                        help="Column separator of --tpm_file, e.g. ',' or 'tab'. Default: ','.")
    parser.add_argument("--biomart_file", type=str, required=True,
                        help="Biomart table with 'Gene stable ID' and 'Transcript stable ID' (data/biomart.txt).")
    parser.add_argument("--regulator_file", type=str, required=True, help="Regulator list (data/regulator_db.tsv).")
    parser.add_argument("--regulator_types", type=str, nargs='+', default=['TF', 'SF', 'TF_SF'],
                        help="Regulator_type values whose genes are regulators.")
    parser.add_argument("--output", type=str, required=True, help="Output network TSV.")
    parser.add_argument("--min_tpm", type=float, default=1.0, help="Minimal TPM of an expressed transcript.")
    parser.add_argument("--min_fraction", type=float, default=0.25,
                        help="Minimal fraction of samples in which a transcript is expressed.")
    parser.add_argument("--benchmark", type=str, required=False,
                        help="Optional TSV with wall time, CPU time and peak memory (this process and its "
                             "dask workers) of each step.")
    parser.add_argument("--num_permutations", type=int, default=1000)
    parser.add_argument("--num_target_clusters", type=int, default=1)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    benchmark = Benchmark()
    with benchmark.step('start_cluster') as info:
        client = Client(LocalCluster())
        info.update(n_workers=len(client.scheduler_info()['workers']), n_threads=sum(client.nthreads().values()))
    network = gene_count_grn(args.tpm_file, args.biomart_file, args.regulator_file, args.regulator_types,
                             args.min_tpm, args.min_fraction, args.num_permutations, args.num_target_clusters,
                             args.seed, client, benchmark,
                             tpm_sep='\t' if args.tpm_sep in ('tab', '\\t') else args.tpm_sep)
    client.close()

    with benchmark.step('write_output'):
        network.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(network)} edges to {args.output}")
    if args.benchmark is not None:
        benchmark.table().to_csv(args.benchmark, sep='\t', index=False)


if __name__ == "__main__":
    main()
