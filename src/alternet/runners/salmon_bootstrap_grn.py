import argparse

import pandas as pd
from distributed import Client, LocalCluster
from sklearn.preprocessing import StandardScaler

from alternet.benchmark import Benchmark
from alternet.inference import compute_grn
from alternet.informative_exons import is_expressed
from alternet.salmon_bootstraps import read_salmon_bootstraps, read_salmon_tpm, salmon_sample_name


def salmon_bootstrap_grn(quant_dirs, regulator_file, regulator_types=('TF', 'SF', 'TF_SF'), min_tpm=1.0,
                         min_fraction=0.25, num_permutations=1000, num_target_clusters=1, seed=11, client=None,
                         benchmark=None, max_fdr_features=20000):
    """
    Transcript GRN from Salmon TPM with signifikante P-values and bootstrap SHAP statistics from the
    Salmon bootstrap replicates (salmon quant --numBootstraps).

    Parameters:
        quant_dirs (list): Salmon output directories (<sample>/ with quant.sf and aux_info/), one per
            sample. The directory name is the sample name.
        regulator_file (str): Regulator list (data/regulator_db.tsv).
        regulator_types (list): Regulator_type values whose transcripts are regulators.
        min_tpm (float): Minimal TPM of an expressed transcript.
        min_fraction (float): Minimal fraction of samples in which a transcript is expressed.
        num_permutations (int): Number of permutations for the empirical P-values.
        num_target_clusters (int): Number of target clusters for the permutation background.
        seed (int): Random seed for the regressors.
        client (dask.distributed.Client): Dask client. A new local cluster is started if None.
        benchmark (Benchmark): Records run time and memory of each step. A new one is used if None.
        max_fdr_features (int): With more expressed transcripts than this, only grnboost2 with the
            bootstraps is run, without the permutations (no 'pvalue' / 'pvalue_bh' columns).

    Returns:
        pd.DataFrame: Network with the columns of compute_grn, including the bootstrap columns.
    """
    if benchmark is None:
        benchmark = Benchmark()

    with benchmark.step('read_inputs') as info:
        quant_dirs = {salmon_sample_name(d): d for d in quant_dirs}
        tpm = read_salmon_tpm(quant_dirs)
        expressed = list(tpm.index[is_expressed(tpm, min_tpm, min_fraction)])
        tpm = tpm.loc[expressed].T
        bootstraps = read_salmon_bootstraps(quant_dirs, expressed)
        print(f"{tpm.shape[0]} samples, {tpm.shape[1]} expressed transcripts, {len(bootstraps)} bootstraps")

        # The models are trained on the standardized point estimates, so the bootstraps get the same scaling
        scaler = StandardScaler().fit(tpm)
        data = pd.DataFrame(scaler.transform(tpm), index=tpm.index, columns=tpm.columns)
        bootstraps = [pd.DataFrame(scaler.transform(b.loc[tpm.index, tpm.columns]), index=tpm.index, columns=tpm.columns)
                      for b in bootstraps]

        regulators = pd.read_csv(regulator_file, sep='\t')
        regulator_transcripts = list(set(regulators.loc[regulators['Regulator_type'].isin(regulator_types),
                                                        'Transcript stable ID']) & set(data.columns))
        print(f"{len(regulator_transcripts)} expressed regulator transcripts")
        info.update(n_samples=data.shape[0], n_transcripts=data.shape[1], n_bootstraps=len(bootstraps))

    run_fdr = data.shape[1] <= max_fdr_features
    if not run_fdr:
        print(f"{data.shape[1]} expressed transcripts > {max_fdr_features}: running only the bootstraps, no FDR")
    with benchmark.step('inference', n_regulators=len(regulator_transcripts), n_targets=data.shape[1],
                        num_permutations=num_permutations if run_fdr else 0, n_bootstraps=len(bootstraps),
                        run_fdr=run_fdr) as info:
        network = compute_grn(data, 'all', regulator_transcripts, client=client, seed=seed,
                              num_target_cluster=num_target_clusters, num_permutations=num_permutations,
                              bootstraps=bootstraps, run_fdr=run_fdr)
        info.update(n_edges=len(network))
    return network


def main():
    parser = argparse.ArgumentParser(
        description="Transcript GRN from Salmon TPM with signifikante P-values and bootstrap SHAP statistics "
                    "from the Salmon bootstrap replicates (salmon quant --numBootstraps)."
    )
    parser.add_argument("--quant_dirs", type=str, nargs='+', required=True,
                        help="Salmon output directories (<sample>/ with quant.sf and aux_info/), one per sample. "
                             "The directory name is the sample name.")
    parser.add_argument("--regulator_file", type=str, required=True, help="Regulator list (data/regulator_db.tsv).")
    parser.add_argument("--regulator_types", type=str, nargs='+', default=['TF', 'SF', 'TF_SF'],
                        help="Regulator_type values whose transcripts are regulators.")
    parser.add_argument("--output", type=str, required=True, help="Output network TSV.")
    parser.add_argument("--min_tpm", type=float, default=1.0, help="Minimal TPM of an expressed transcript.")
    parser.add_argument("--min_fraction", type=float, default=0.25,
                        help="Minimal fraction of samples in which a transcript is expressed.")
    parser.add_argument("--benchmark", type=str, required=False,
                        help="Optional TSV with wall time, CPU time and peak memory (this process and its "
                             "dask workers) of each step.")
    parser.add_argument("--max_fdr_features", type=int, default=20000,
                        help="With more expressed transcripts than this, only the bootstraps are run, without the "
                             "permutations (no pvalue / pvalue_bh columns). Default: 20000.")
    parser.add_argument("--num_permutations", type=int, default=1000)
    parser.add_argument("--num_target_clusters", type=int, default=1)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    benchmark = Benchmark()
    with benchmark.step('start_cluster') as info:
        client = Client(LocalCluster())
        info.update(n_workers=len(client.scheduler_info()['workers']), n_threads=sum(client.nthreads().values()))
    network = salmon_bootstrap_grn(args.quant_dirs, args.regulator_file, args.regulator_types, args.min_tpm,
                                   args.min_fraction, args.num_permutations, args.num_target_clusters, args.seed,
                                   client, benchmark, args.max_fdr_features)
    client.close()

    with benchmark.step('write_output'):
        network.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(network)} edges to {args.output}")
    if args.benchmark is not None:
        benchmark.table().to_csv(args.benchmark, sep='\t', index=False)


if __name__ == "__main__":
    main()
