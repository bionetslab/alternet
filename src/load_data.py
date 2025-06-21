from tqdm import tqdm
import pandas as pd
import re



def read_tf_list(tf_path, biomart):
    ''' 
    Reads transcription factor list from file and merges it with gene information from biomart dataframe

    Attributes
    ---------------
    tf_path: string, path to tab seperated file of transcription factors
    biomart: Dataframe

    Returns
    ---------------
    tf_list: Dataframe which containes transcription factors with their corresponding gene and trancsript information from biomart
    '''

    tf_list = pd.read_csv(tf_path, sep='\t', header = None)
    tf_list.columns  = ['TF'] # rename column to TF 
    tf_list = tf_list.merge(biomart, left_on = 'TF', right_on = 'Gene name') # merge tf list with biomart: join 'TF' with 'Gene name'
    tf_list = tf_list.loc[:, ['TF', 'Gene stable ID', 'Transcript stable ID']].drop_duplicates() #remove individual versions
    tf_list = tf_list.reset_index()
    return tf_list


class TissueNotFoundException(Exception):
    '''
    custom exception
    used when requested Tissue type is not found in GTEX data
    '''
    pass


def retrieve_GTEX_tissue_sampleids(annotation_file, tissue='Liver'):
    ''' 
    Reads the annotation file to filter out samples that belong to the specified tissue.
    If tissue is not found it raises an TissueNotFoundException

    Attributes
    ---------------
    annotation_file: str; path to annotation file
    tissue: str, tissue type


    Returns
    ---------------
    A list of sample IDs corresponding to the selected tissue.
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
    Reads a subset of columns from the GTEX transcript expression file based on sample IDs specific to tissue type.

    Attributes
    ---------------
    path: str; Path to the GTEX transcript expression data
    sample_ids: List; List of sample IDs for specified tissue
    headers: List, List of mandatory columns to read

    Returns
    ---------------
    Dataframe; containing the specified transcript or gene expression data for the given sample IDs
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
    Reads a subset of columns from the GTEX transcript expression file based on sample IDs specific to tissue type.

    Attributes
    ---------------
    path: str; Path to the GTEX transcript expression data
    sample_ids: List; List of sample IDs for specified tissue
    gene_ids: gene_ids used for grn
    headers: List, List of mandatory columns to read

    Returns
    ---------------
    Dataframe; containing the specified transcript or gene expression data for the given sample IDs
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
    Removes version numbers of transcript IDs 

    Attributes
    ---------------
    data: Dataframe; contains ID column
    transcript_column: str, name of column with transcript IDs

    Returns
    ---------------
    Dataframe; with version number removed from transcript IDs
    '''   
    data[transcript_column] = data[transcript_column].str.split('.').str[0]

    return data

def clean_GTEX_tissue_transcript_counts(data, biomart, relevant_columns=['transcript_id', 'gene_id'], biomart_column='Transcript stable ID'):

                    
    ''' 
    pipeline for removing version IDs of Transcripts and  Genes,   filtering for protein-coding genes, and removing genes with low expression


    Attributes
    ---------------
    data: Dataframe; transcript expression data
    biomart: Dataframe contains transcript annotations from Biomart
    biomart_column: Column name in biomart to use for filtering
    relevant_columns: List, List of columns to clean
    biomart_columns: List, List of columns to use for filtering

    Returns
    ---------------
    Dataframe; transcript expression data, filtered for protein-coding genes with gene_id as the index
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



def separate_tf_genes(data, tf_list, data_column='transcript_id', biomart_column='Transcript stable ID'):
    ''' 
    Separates gene expression data for transcription factors and genes
    
    Attributes
    ---------------
    data: Dataframe; gene expression data
    tf_list: Dataframe; contains transcription factor information
    data_column: str, name of column in gene expression data
    biomart_column: str, name of column in tf_list

    Returns
    ---------------
    Dataframe; containing gene expression data for transcription factors
    '''   
    tfs = data[data[data_column].isin(tf_list[biomart_column])]
    genes = data[~data[data_column].isin(tf_list[biomart_column])]
    assert tfs.shape[0] + genes.shape[0] == data.shape[0]

    return tfs, genes



def get_target_genes(data):
    ''' 
    Gets target genes from gene expression data

    Attributes
    ---------------
    data: Dataframe; gene expression data

    Returns
    ---------------
    List; containing gene IDs
    '''   
    return data.index.tolist()


def load_data(config, biomart, tf_list):
    '''
    Loads gene expression data and transcription factor list
    
    Attributes
    --------------- 
    config: dict; configuration dictionary containing paths to data files

    Returns
    ---------------
    transcript_tfs: Dataframe; containing transcription factor transcript expression data of tfs
    gene_tfs: Dataframe; containing transcription factor gene expression data of tfs
    targets: Dataframe; containing gene expression data of genes that are not tfs

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


def prepare_for_inference(transcript_tfs, gene_tfs, targets, transcript_column='transcript_id', gene_column='gene_id'):
    '''
    Prepares gene expression data for inference by separating it into canonical and as-aware gene expression data
    
    Attributes
    ---------------
    tfs: Dataframe; contains transcription factor gene expression data
    genes: Dataframe; contains gene expression data


    Returns
    ---------------
    data_canonical: Dataframe; containing canonical gene expression data    
    data_asaware: Dataframe; containing as-aware gene expression data
    target_gene_list: List; containing target gene IDs
    
    '''

    print('Preparing data for inference: Separation into canonical and as-aware gene expression data')


    
    data_canonical = pd.concat([gene_tfs, targets], axis=0)
    assert data_canonical.shape[0] == targets.shape[0] + gene_tfs.shape[0]
    assert data_canonical.shape[1] == targets.shape[1]
    data_canonical = data_canonical.set_index(gene_column)

    # get as-aware gene expression data
    tf_isoform = transcript_tfs.copy(deep=True)
    tf_isoform[gene_column] = tf_isoform[transcript_column]
    tf_isoform = tf_isoform.drop(columns=transcript_column)
    tf_isoform = tf_isoform.set_index(gene_column)

    data_asware = pd.concat([tf_isoform, targets.set_index(gene_column)], axis=0)
    assert data_asware.shape[0] == targets.shape[0] + tf_isoform.shape[0]
    assert data_asware.shape[1] + 1 == targets.shape[1] #+1 because setting index removes column


    target_gene_list = targets[gene_column].unique().tolist()

    return data_canonical, data_asware, target_gene_list


