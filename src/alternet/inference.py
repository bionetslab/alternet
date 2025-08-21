
import pandas as pd
import os.path as op
import os
from distributed import Client, LocalCluster
from arboreto_added.algo import grnboost2
from data_preprocessing import *
from collections import defaultdict
from tqdm import tqdm 


def compute_grn(gene_data, target_names, tf_list, client, use_tf=True):
    
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
    
    print('Computing network')

    # compute the GRN
    if not use_tf:
        network = grnboost2(expression_data=gene_data, 
                            client_or_address=client)
    else:
        network = grnboost2(expression_data=gene_data,
                            target_genes=target_names,
                            tf_names=tf_list,
                            client_or_address=client)

    # write the GRN to file
    network.columns = ['source', 'target', 'importance']
    return network



def inference(config, nruns, data_canonical, data_asware, target_gene_list, tf_list, aggregate=True):
    '''
    Performs inference to create gene regulatory networks (GRNs) for transcript-level and gene-level data.
    Optionally aggregates the results from multiple runs.

    Parameters:
        config (dict): Configuration dictionary containing paths and settings for the inference process.
        nruns (int): Number of inference runs to perform.
        aggregate (bool): Whether to aggregate results from multiple runs. Default is True.

    Returns:
        tuple:
            - as_aware_grn (pd.DataFrame): Inferred or aggregated AS-aware GRN.
            - canonical_grn (pd.DataFrame): Inferred or aggregated canonical GRN.
    
    '''

    print('Starting inference ...')

    client = Client(LocalCluster())


    
    ## RUN INFERENCE for AS-Aware network
    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')
    os.makedirs(results_dir_grn, exist_ok=True)

    data_asaware_t = data_asware.T
    data_canonical_t = data_canonical.T

    tfs_transcripts = tf_list['Transcript stable ID'].unique().tolist()
    tfs_genes = tf_list['Gene stable ID'].unique().tolist()
    
    as_aware_grns = []
    for i in range(1, nruns+1):
        grn1 = compute_grn(data_asaware_t,
                                    target_gene_list,
                                    tfs_transcripts,
                                    client,
                                    use_tf=True)
        as_aware_grns.append(grn1)

    canoncial_grns = []
    for i in range(1, nruns+1):
        grn2 = compute_grn(data_canonical_t,
                                    target_gene_list,
                                    tfs_genes,
                                    client,
                                    use_tf=True)
        canoncial_grns.append(grn2)
        
    
    client.close()
    print('Inference complete')


    if aggregate:
        print('Aggregate results')

        as_aware = aggregate_results(as_aware_grns)
        canoncial = aggregate_results(canoncial_grns)
        return as_aware, canoncial
    
    return as_aware_grns, canoncial_grns




def aggregate_results(grn_results):
    '''
    Aggregates results from multiple GRN inference runs.

    Parameters:
        grn_results (list of pd.DataFrame): Results from multiple GRNBoost runs.

    Returns:
        pd.DataFrame: Aggregated consensus network.
    '''
    # Create a dictionary to store the aggregated results
    edge_data = defaultdict(list)

    print("Going through edges of each df ...")
    for df in grn_results:
        for _, row in tqdm(df.iterrows(), total=df.shape[0]):
            edge = (row['source'], row['target'])
            edge_data[edge].append(row['importance'])

    aggregated_edges = []
    for (source, target), weights in tqdm(edge_data.items(), total=len(edge_data)):
        aggregated_edges.append(
            {
                'source': source,
                'target': target,
                'frequency': len(weights),
                'mean_importance' : sum(weights) / len(weights),
                "median_importance": sorted(weights)[len(weights) // 2]
            }
        )
    aggregated_df = pd.DataFrame(aggregated_edges)

    return aggregated_df

    