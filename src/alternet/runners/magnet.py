import pandas as pd
from inference import *
from load_data import *
from total_pipeline import *

def load_magnet_data(config, tf_list):
    '''
    Loads MAGNet gene and transcript expression data filtered for a specific tissue and separates transcription factors from target genes.

    Parameters:
        config (dict): Configuration dictionary containing paths to data files and tissue name. Expected keys:
            - 'transcript_data' (str): Path to transcript-level expression data.
            - 'count_data' (str): Path to gene-level expression data.
            - 'sample_attributes' (str): Path to sample attribute file.
            - 'tissue' (str): Tissue name to filter samples.
        tf_list (pd.DataFrame): DataFrame containing transcription factor information.

    Returns:
        tuple:
            - transcript_tfs (pd.DataFrame): Transcript-level expression data for transcription factors.
            - gene_tfs (pd.DataFrame): Gene-level expression data for transcription factors.
            - targets (pd.DataFrame): Gene-level expression data for target genes (non-transcription factors).

    
    '''
    #transcripts already preprocessed
    transcript_data = pd.read_csv(config['transcript_data'], index_col=0)
    gene_ids = transcript_data['gene_id'].unique()
    gene_data = pd.read_csv(config['count_data'], index_col=0)
    gene_data.index.name = 'gene_id'
    gene_data = gene_data.reset_index()

    gene_data = gene_data[gene_data['gene_id'].isin(gene_ids)]

    sample_attributes = pd.read_csv(config['sample_attributes'])
    sample_attributes = sample_attributes.loc[:, ['sample_name', 'etiology']]
    samples = sample_attributes[sample_attributes['etiology'] == config['tissue']]
    samples = samples['sample_name'].tolist()

    gene_data = gene_data.loc[:, ['gene_id'] + samples ]
    transcript_data = transcript_data.loc[:,['gene_id', 'transcript_id'] + samples]

    #### Inference for HCM
    gene_tfs, targets = separate_tf_genes(gene_data, tf_list, data_column='gene_id', biomart_column='Gene stable ID')
    transcript_tfs, __ = separate_tf_genes(transcript_data, tf_list)

    return transcript_tfs, gene_tfs, targets
    

def main():
    config_file = "/data/bionets/mi34qaba/SpliceAwareGRN/configs/MAGNet_NF.yaml"


    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    # important files to load
    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)


    transcript_tfs, gene_tfs, targets = load_magnet_data(config, tf_list)

    inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets)

    
if __name__ == "__main__":
    main()