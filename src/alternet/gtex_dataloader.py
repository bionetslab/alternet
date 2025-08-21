import pandas as pd
from data_preprocessing import *

class TissueNotFoundException(Exception):
    '''
    custom exception
    used when requested Tissue type is not found in GTEX data
    '''
    pass


def retrieve_GTEX_tissue_sampleids(annotation_file, tissue='Liver'):
    ''' 
    Reads the annotation file to filter out samples that belong to the specified tissue.
    Raises a TissueNotFoundException if the tissue is not found.

    Parameters:
        annotation_file (str): Path to the annotation file.
        tissue (str): Tissue type to filter for.

    Returns:
        list: Sample IDs corresponding to the selected tissue.
    '''

    print('Retrieving tissue sample IDs')
    sample_meta = pd.read_csv(annotation_file, sep='\t')
    sub = sample_meta.loc[:, ['SAMPID', 'SMTS']] # SMTS == tissue type
    try:
        sub = sub[sub.SMTS == tissue]
    except:
        raise TissueNotFoundException(f'{tissue} is not a valid tissue identifier in annotation file')
    return sub.SAMPID.tolist()


def read_GTEX_transcript_expression(path, sample_ids, headers=['transcript_id', 'gene_id']):
    ''' 
    Reads a subset of columns from the GTEx transcript expression file based on sample IDs specific to a tissue type.

    Parameters:
        path (str): Path to the GTEx transcript expression data.
        sample_ids (list): List of sample IDs for the specified tissue.
        headers (list): List of mandatory columns to read.

    Returns:
        pd.DataFrame: Expression data for the specified transcripts or genes and sample IDs.
    '''    
    print('Reading Transcript expression data')
    columns_to_read = headers + sample_ids 
    data = pd.read_csv(path, sep='\t', comment='#', skiprows=2, skipinitialspace=True, nrows=1, header=0)
    columns = data.columns
    columns_to_read = set(columns_to_read).intersection(set(columns))
    data = pd.read_csv(path, sep='\t', comment='#', skiprows=2, skipinitialspace=True, usecols=columns_to_read,header=0)   

    return data

def read_GTEX_gene_expression(path, sample_ids, gene_ids, headers=['Name']):
    ''' 
    Reads a subset of columns from the GTEx transcript expression file based on sample IDs specific to a tissue type.

    Parameters:
        path (str): Path to the GTEx transcript expression data.
        sample_ids (list): List of sample IDs for the specified tissue.
        gene_ids (list): Gene IDs used for GRN inference.
        headers (list): List of mandatory columns to read.

    Returns:
        pd.DataFrame: Expression data for the specified transcripts or genes and sample IDs.
    '''    
    print('Reading Gene expression data')
    columns_to_read = headers + sample_ids 
    data = pd.read_csv(path, sep='\t', comment='#', skiprows=2, skipinitialspace=True, nrows=1, header=0)
    columns = data.columns
    columns_to_read = set(columns_to_read).intersection(set(columns))
    data = pd.read_csv(path, sep='\t', comment='#', skiprows=2, skipinitialspace=True, usecols=columns_to_read,header=0)
    
    data = remove_version_id(data, transcript_column='Name')
    data.rename(columns={"Name": "gene_id"}, inplace=True)

    data = data[data['gene_id'].isin(gene_ids)]

    return data



def remove_version_id(data, transcript_column='transcript_id'):
    ''' 
    Removes version numbers from transcript IDs.

    Parameters:
        data (pd.DataFrame): DataFrame containing a column with transcript IDs.
        transcript_column (str): Name of the column containing transcript IDs.

    Returns:
        pd.DataFrame: DataFrame with version numbers removed from transcript IDs.
    '''   
    data[transcript_column] = data[transcript_column].str.split('.').str[0]

    return data

def clean_GTEX_tissue_transcript_counts(data, biomart, relevant_columns=['transcript_id', 'gene_id'], biomart_column='Transcript stable ID'):

                    
    ''' 
    Pipeline for removing version numbers from transcript and gene IDs, filtering for protein-coding genes,
    and removing genes with low expression.

    Parameters:
        data (pd.DataFrame): Transcript expression data.
        biomart (pd.DataFrame): DataFrame containing transcript annotations from BioMart.
        biomart_column (str): Column name in BioMart used for filtering.
        relevant_columns (list): List of columns to clean (e.g., remove version numbers).
        biomart_columns (list): List of columns to use for filtering.

    Returns:
        pd.DataFrame: Filtered transcript expression data with gene_id as the index.
    '''   

    
    print('Cleaning up counts')

    for i in range(len(relevant_columns)):
        data = remove_version_id(data, transcript_column=relevant_columns[i])

    # filter for protein coding genes in biomart and select expression values of these in gex
    data = data[
        data[relevant_columns[0]].isin(biomart[biomart['Gene type'] == 'protein_coding'][biomart_column].tolist())]
    
    threshold = data.shape[1] * 0.1
    mask = (data == 0).sum(axis=1) > threshold
    data = data[~mask]
    
    return data

def load_gtex_data(config, biomart, tf_list):
    '''
    Loads gene expression data and transcription factor list.

    Parameters:
        config (dict): Configuration dictionary containing paths to data files.

    Returns:
    
        - transcript_tfs (pd.DataFrame): Transcription factor transcript expression data.
        - gene_tfs (pd.DataFrame): Transcription factor gene expression data.
        - targets (pd.DataFrame): Gene expression data of genes that are not transcription factors.

    '''

    # Retrieve tissue sample IDs
    
    tissue_ids = retrieve_GTEX_tissue_sampleids(config['sample_attributes'], tissue=config['tissue'])
    

    # load transcript data
    transcript_data = read_GTEX_transcript_expression(config['transcript_data'], tissue_ids)
    transcript_data = clean_GTEX_tissue_transcript_counts(transcript_data, biomart)


    # load gene data
    gene_ids = transcript_data['gene_id'].unique()
    gene_data = read_GTEX_gene_expression(config['count_data'], tissue_ids, gene_ids, headers=['Name'])


    gene_tfs, targets = separate_tf_genes(gene_data, tf_list, data_column='gene_id', biomart_column='Gene stable ID')


    # separate transcription factors and genes
    transcript_tfs, __ = separate_tf_genes(transcript_data, tf_list)


    return transcript_tfs, gene_tfs, targets
