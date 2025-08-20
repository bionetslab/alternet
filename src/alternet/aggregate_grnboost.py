import pandas as pd
from collections import defaultdict
import argparse
from tqdm import tqdm 
import os.path as op


def read_grn_files(path, file_prefix, nruns):
    '''
    Loads all GRN repetitions of a given type.

    Parameters:
        config (dict): Configuration dictionary containing:
            - 'results_dir' (str): Base directory where results are stored.
            - 'tissue' (str): Tissue name used as a prefix to locate GRN result files.
            - 'nruns': number of repetitions
    Returns:
        list: List of all loaded GRN repetitions.
    '''
    grn_results = []
    for i in tqdm(range(1, nruns+1)):
        data = pd.read_csv(op.join(path, file_prefix + str(i) + '.tsv'), sep='\t')
        grn_results.append(data)

    return grn_results

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

    