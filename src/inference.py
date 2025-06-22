
import pandas as pd
import os.path as op
import os
from distributed import Client, LocalCluster
from arboreto_added.algo import grnboost2
import yaml
import argparse
from load_data import *
from aggregate_grnboost import *


def compute_and_save_network(gene_data, target_names, tf_list, client, file, use_tf=False, aggregate=True):
    
    ''' 
    computes GRN using GRNBoost2 
    
    Attributes
    ---------------
    data: Dataframe; expression data 
    tf_list: List, List of transcription factors to be included
    client: Dask client; used to distribute computation across a cluster
    file: str, output file path where GRN will be saved
    use_tf: boolean flag; if true, list of transcription factors is used during network computation. If False, all genes are used.

    Returns
    ---------------
    The function returns the computed network, a DataFrame that contains the regulatory interactions between genes and/or transcription factors.
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
    network.to_csv(file, sep='\t', index=False, header=True)
    
    return network
    

def aggregate_and_save(results_dir_grn, prefix, nruns, type='as-aware_'):
    grn = read_grn_files(results_dir_grn, prefix, nruns)
    aggregate_df = aggregate_results(grn)
    aggregate_df.to_csv(op.join(results_dir_grn, prefix + 'aggregated.tsv'), sep='\t', index=False)
    return aggregate_df


def inference(config, nruns, data_canonical, data_asware, target_gene_list, tf_list, aggregate=True):
    """
    Perform inference to create gene regulatory networks (GRNs) for transcript data and gene data.
    Optionally, it can aggregate the results from multiple runs.

    Attributes:
    config (dict): Configuration dictionary containing paths and settings for the inference process.
    nruns (int): Number of runs to perform for the inference.
    aggregate (bool): Whether to aggregate the results from multiple runs. Default is True.
    
    Returns:
    None
    """



    print('Starting inference ...')


    client = Client(LocalCluster())


    
    ## RUN INFERENCE for AS-Aware network
    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')
    os.makedirs(results_dir_grn, exist_ok=True)

    as_aware_prefix = f"{config['tissue']}_as-aware.network_"
    canonical_prefix = f"{config['tissue']}_canonical.network_"

    data_asaware_t = data_asware.T
    data_canonical_t = data_canonical.T

    tfs_transcripts = tf_list['Transcript stable ID'].unique().tolist()
    tfs_genes = tf_list['Gene stable ID'].unique().tolist()
    
    for i in range(1, nruns+1):
        print(i)
        print()
        file = op.join(results_dir_grn, as_aware_prefix + f"{i}.tsv")
        grn1 = compute_and_save_network(data_asaware_t,
                                    target_gene_list,
                                    tfs_transcripts,
                                    client,
                                    file,
                                    use_tf=True)
        

    for i in range(1, nruns+1):
        print(i)
        print()
        file = op.join(results_dir_grn, canonical_prefix + f"{i}.tsv")
        grn2 = compute_and_save_network(data_canonical_t,
                                    target_gene_list,
                                    tfs_genes,
                                    client,
                                    file,
                                    use_tf=True)
        
    
    client.close()
    print('Inference complete')
    print('Results saved to:', results_dir_grn)


    if aggregate:
        print('Aggregate results')

        as_aware = aggregate_and_save(results_dir_grn, as_aware_prefix, nruns, type='as-aware_')
        canoncial = aggregate_and_save(results_dir_grn, canonical_prefix, nruns, type='canonical_')
        print('Aggregated results saved to:', results_dir_grn)
        return as_aware, canoncial
    
    return grn1, grn2

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process a file from the command line.")
    
    # Add the file argument
    parser.add_argument('-f', type=str, help='The file to process')
    parser.add_argument('-n', type=int, help='Number of runs of grnboost')

    # Parse the arguments
    args = parser.parse_args()
    with open(args.f, 'r') as f:
        config = yaml.safe_load(f)

    
    print('Load data and prepare for inference')
    
    # important files to load
    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)
    
    
    transcript_tfs, gene_tfs, targets = load_data(config, biomart, tf_list)
    data_canonical, data_asware, target_gene_list = prepare_for_inference(transcript_tfs, gene_tfs, targets)


    inference(config, args.n, data_canonical=data_canonical, data_asware=data_asware, target_gene_list=target_gene_list, tf_list=tf_list)