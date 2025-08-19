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

def main():
    parser = argparse.ArgumentParser(description="Process a file from the command line.")
    # Add the file argument
    parser.add_argument('-f', type=str, help='The file to process')
    parser.add_argument('-p', type=str, help='Prefix of the files')
    parser.add_argument('-n', type=int, help='Number of runs of grnboost')
    args = parser.parse_args()
    # Read GRNBoost files
    print("Reading files ...")
    grn_results = read_grn_files(args.f, args.p, args.n)
    # Aggregate results
    print("Aggregating results ...")
    aggregated_df = aggregate_results(grn_results)
    # Save aggregated results
    print("Saving results ...")
    aggregated_df.to_csv(op.join(args.f, args.p + 'aggregated.tsv'), sep='\t', index=False)
    print('Saved results in ', op.join(args.f, args.p + 'aggregated.tsv'))
    print("Done!")
if __name__ == "__main__":
    main()
    