import pandas as pd


def separate_tf_genes(data, tf_list, data_column='transcript_id', biomart_column='Transcript stable ID'):
    ''' 
    Separates gene expression data into transcription factors and other genes.

    Parameters:
        data (pd.DataFrame): Gene expression data.
        tf_list (pd.DataFrame): DataFrame containing transcription factor information.
        data_column (str): Name of the column in the gene expression data to match.
        biomart_column (str): Name of the column in tf_list to match.

    Returns:
        pd.DataFrame: Gene expression data for transcription factors.
    '''   
    tfs = data[data[data_column].isin(tf_list[biomart_column])]
    genes = data[~data[data_column].isin(tf_list[biomart_column])]
    assert tfs.shape[0] + genes.shape[0] == data.shape[0]

    return tfs, genes



def get_target_genes(data):
    ''' 
    Gets target genes from gene expression data.

    Parameters:
        data (pd.DataFrame): Gene expression data.

    Returns:
        list: List containing gene IDs.
    '''   
    return data.index.tolist()




def prepare_for_inference(transcript_tfs, gene_tfs, targets, transcript_column='transcript_id', gene_column='gene_id'):
    '''
    Prepares gene expression data for inference by separating it into canonical and AS-aware gene expression data.

    Parameters:
        transcript_tfs (pd.DataFrame): Transcription factor isoform expression data
        gene_tfs (pd.DataFrame): Transcription factor gene expression data.
        targets (pd.DataFrame): Gene expression data of non-tf gene expression data.
        transcript_column (str): column header of transcript ID column
        gene_column (str): column header of gene ID column

    Returns:
        tuple:
            - data_canonical (pd.DataFrame): Canonical gene expression data.
            - data_asaware (pd.DataFrame): AS-aware gene expression data.
            - target_gene_list (list): List of target gene IDs.
    
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


