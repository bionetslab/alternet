import pandas as pd
import os.path as op
import os
from distributed import Client, LocalCluster
from signifikante.algo import signifikante_fdr, grnboost2
from alternet.data_preprocessing import *
from collections import defaultdict
from tqdm import tqdm


def compute_grn(gene_data, target_names, tf_list, client=None, seed = 11, num_target_cluster = 1, num_permutations=1000,
                bootstraps=None, apply_bh_correction=True, run_fdr=True):

    '''
    Computes a gene regulatory network (GRN) with GRNBoost2 and empirical P-values (signifikante_fdr).

    Parameters:
        gene_data (pd.DataFrame): Expression data, samples x features. All columns are used as targets
            during inference.
        target_names (list or 'all'): Targets to keep in the returned network.
        tf_list (list): Regulators, a subset of the columns of gene_data.
        client (dask.distributed.Client): Dask client. A new local cluster is started if None.
        seed (int): Random seed for the regressors.
        num_target_cluster (int): Number of target clusters for the permutation background.
        num_permutations (int): Number of permutations for the empirical P-values.
        bootstraps (list): Optional bootstrapped expression matrices with the same shape and columns as
            gene_data. Adds the columns 'shap_mean', 'shap_median', 'shap_std' and 'max_rank_diff'.
        apply_bh_correction (bool): Adds Benjamini-Hochberg adjusted P-values in 'pvalue_bh'.
        run_fdr (bool): If False, only runs grnboost2 (with the bootstraps), without the permutations.
            The network then has no 'pvalue' / 'pvalue_bh' columns.

    Returns:
        pd.DataFrame: Network with columns 'source', 'target', 'importance', 'pvalue' (and 'pvalue_bh',
            bootstrap columns if requested).
    '''

    if client is None:
        client = Client(LocalCluster())
    if not run_fdr:
        # the same grnboost2 call that signifikante_fdr makes before its permutations
        network = grnboost2(expression_data=gene_data, tf_names=tf_list, client_or_address=client, seed=seed,
                            bootstraps=bootstraps)
    else:
        # signifikante infers the GRN for all columns of gene_data and adds empirical p-values
        network = signifikante_fdr(
            expression_data=gene_data,
            cluster_representative_mode='medoid',
            num_target_clusters=num_target_cluster,
            target_cluster_mode='kmeans',
            tf_names=tf_list,
            client_or_address=client,
            seed=seed,
            num_permutations=num_permutations,
            apply_bh_correction=apply_bh_correction,
            bootstraps=bootstraps
        )

    network = network.rename(columns={'TF': 'source'})

    # signifikante has no target subset for clustered modes, so restrict targets afterwards
    if target_names != 'all':
        network = network[network['target'].isin(set(target_names))]

    return network.reset_index(drop=True)
