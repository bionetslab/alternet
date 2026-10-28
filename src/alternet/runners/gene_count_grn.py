import argparse

import pandas as pd
from distributed import Client, LocalCluster

from alternet.data_preprocessing import standardize_dataframe
from alternet.inference import compute_grn
from alternet.informative_exons import is_expressed


def read_gene_counts(counts_file):
    """
    Reads a merged gene count table (nf-core/rnaseq salmon.merged.gene_counts.tsv: gene_id, gene_name,
    one column per sample). Non-numeric columns are dropped and the version suffix is removed from the
    gene IDs; genes that share an ID without version (e.g. _PAR_Y copies) are summed.

    Returns:
        pd.DataFrame: gene_ids (index) x samples.
    """
    counts = pd.read_csv(counts_file, sep='\t', index_col='gene_id').select_dtypes('number')
    counts.index = counts.index.str.split('.').str[0]
    return counts.groupby(level=0).sum()


def gene_count_grn(counts_file, regulator_file, regulator_types=('TF', 'SF', 'TF_SF'), min_cpm=1.0,
                   min_fraction=0.25, num_permutations=1000, num_target_clusters=1, seed=11, client=None):
    """
    Gene GRN from gene counts (normalized to CPM) with regulator genes as regulators, with
    signifikante P-values.

    Parameters:
        counts_file (str): Merged gene count table (see read_gene_counts).
        regulator_file (str): Regulator list (data/regulator_db.tsv).
        regulator_types (list): Regulator_type values whose genes are regulators.
        min_cpm (float): Minimal counts per million of an expressed gene.
        min_fraction (float): Minimal fraction of samples in which a gene is expressed.
        num_permutations (int): Number of permutations for the empirical P-values.
        num_target_clusters (int): Number of target clusters for the permutation background.
        seed (int): Random seed for the regressors.
        client (dask.distributed.Client): Dask client. A new local cluster is started if None.

    Returns:
        pd.DataFrame: Network with the columns of compute_grn.
    """
    counts = read_gene_counts(counts_file)
    cpm = counts / counts.sum() * 1e6
    cpm = cpm.loc[is_expressed(cpm, min_cpm, min_fraction)].T
    print(f"{cpm.shape[0]} samples, {cpm.shape[1]} expressed genes")
    data = standardize_dataframe(cpm)

    regulators = pd.read_csv(regulator_file, sep='\t')
    regulator_genes = list(set(regulators.loc[regulators['Regulator_type'].isin(regulator_types),
                                              'Gene stable ID']) & set(data.columns))
    print(f"{len(regulator_genes)} expressed regulator genes")

    return compute_grn(data, 'all', regulator_genes, client=client, seed=seed,
                       num_target_cluster=num_target_clusters, num_permutations=num_permutations)


def main():
    parser = argparse.ArgumentParser(
        description="Gene GRN from gene counts with regulator genes as regulators, with signifikante P-values."
    )
    parser.add_argument("--gene_counts", type=str, required=True,
                        help="Merged gene count table (tab-separated, 'gene_id' column, one column per sample).")
    parser.add_argument("--regulator_file", type=str, required=True, help="Regulator list (data/regulator_db.tsv).")
    parser.add_argument("--regulator_types", type=str, nargs='+', default=['TF', 'SF', 'TF_SF'],
                        help="Regulator_type values whose genes are regulators.")
    parser.add_argument("--output", type=str, required=True, help="Output network TSV.")
    parser.add_argument("--min_cpm", type=float, default=1.0, help="Minimal counts per million of an expressed gene.")
    parser.add_argument("--min_fraction", type=float, default=0.25,
                        help="Minimal fraction of samples in which a gene is expressed.")
    parser.add_argument("--num_permutations", type=int, default=1000)
    parser.add_argument("--num_target_clusters", type=int, default=1)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    client = Client(LocalCluster())
    network = gene_count_grn(args.gene_counts, args.regulator_file, args.regulator_types, args.min_cpm,
                             args.min_fraction, args.num_permutations, args.num_target_clusters, args.seed, client)
    client.close()

    network.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(network)} edges to {args.output}")


if __name__ == "__main__":
    main()
