import pandas as pd


def map_tf_ids(tf_list, biomart):
    ''' 
    Reads a transcription factor list from file and merges it with gene information from a BioMart DataFrame.

    Parameters:
        tf_path (str): Path to a tab-separated file containing transcription factors.
        biomart (pd.DataFrame): DataFrame with gene and transcript information from BioMart.

    Returns:
        pd.DataFrame: Transcription factors with corresponding gene and transcript information from BioMart.
   '''
    tf_list.columns  = ['TF'] # rename column to TF 
    tf_list = tf_list.merge(biomart, left_on = 'TF', right_on = 'Gene name') # merge tf list with biomart: join 'TF' with 'Gene name'
    tf_list = tf_list.loc[:, ['TF', 'Gene stable ID', 'Transcript stable ID']].drop_duplicates() #remove individual versions
    tf_list = tf_list.reset_index()
    return tf_list
