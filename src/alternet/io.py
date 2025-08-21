import pandas as pd
import tqdm 
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


def save_grn_files(grns, path, file_prefix):
    ## Implement