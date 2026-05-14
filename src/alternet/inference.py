
import pandas as pd
import os.path as op
import os
from distributed import Client, LocalCluster
from signifikante.algo import grnboost2
from alternet.data_preprocessing import *
from collections import defaultdict
from tqdm import tqdm 


def compute_grn(gene_data, target_names, tf_list, client=None, use_tf=True, seed = None):
    
    ''' 
    Computes a gene regulatory network (GRN) using GRNBoost2.

    Parameters:
        data (pd.DataFrame): Expression data.
        tf_list (list): List of transcription factors to include.
        client (dask.distributed.Client): Dask client used to distribute computation across a cluster.
        file (str): Output file path where the GRN will be saved.
        use_tf (bool): If True, the transcription factor list is used during network computation. If False, all genes are used.

    Returns:
        pd.DataFrame: Computed network containing regulatory interactions between genes and/or transcription factors.
    '''   

    # compute the GRN
    if not use_tf:
        network = grnboost2(expression_data=gene_data, 
                            client_or_address=client, 
                            seed = seed)
    else:
        network = grnboost2(expression_data=gene_data,
                            target_names=target_names,
                            tf_names=tf_list,
                            client_or_address=client,
                            seed = seed)

    # write the GRN to file
    network.columns = ['source', 'target', 'importance']
    return network




def inference(gene_data,  tf_list, target_names='all', n_runs = 10, set_seed = True):
    '''
    Performs inference to create gene regulatory networks (GRNs) for transcript-level and gene-level data.
    Optionally aggregates the results from multiple runs.

    Parameters:
        config (dict): Configuration dictionary containing paths and settings for the inference process.
        nruns (int): Number of inference runs to perform.
        aggregate (bool): Whether to aggregate results from multiple runs. Default is True.
        set_seed: Automatically sets the seed to the interation index.

    Returns:
        tuple:
            - as_aware_grn (pd.DataFrame): Inferred or aggregated AS-aware GRN.
            - canonical_grn (pd.DataFrame): Inferred or aggregated canonical GRN.
    
    '''


    client = Client(LocalCluster())

    
    grns = []
    for i in tqdm(range(n_runs)):
        if set_seed:
            seed = i
        else:
            seed = None
        grn = compute_grn(gene_data=gene_data,
                            target_names = target_names,
                            tf_list = tf_list,
                            client=client,
                            use_tf=True,
                            seed = seed)
        grns.append(grn)
    
    client.close()

    grn = aggregate_results(grns)

    return grn



def aggregate_results(grn_results):
    '''
    Aggregates results from multiple GRN inference
    Parameters:
        grn_results (list of pd.DataFrame): Results from multiple GRNBoost runs.

    Returns:
        pd.DataFrame: Aggregated consensus network.
    '''
    
    combined_df = pd.concat(grn_results, ignore_index=True)

    aggregated_df = combined_df.groupby(['source', 'target'])['importance'].agg(
        frequency='count',
        mean_importance='mean',
        median_importance='median'
    ).reset_index() 
    
    return aggregated_df



    