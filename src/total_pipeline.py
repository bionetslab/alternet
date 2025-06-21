import pandas as pd
from inference import *
from load_data import *
import utils_database
import utils_network
import json
import yaml
import argparse



def inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets):
    '''
    Executes the full pipeline for GRN inference and annotation with alternative splicing awareness.

    This function performs the following steps:
    1. Loads required resources including biomart annotations and transcription factor list.
    2. Creates a transcript annotation database using APPRIS and Digger data.
    3. Prepares input data for GRN inference based on isoform- and gene-level TFs and targets.
    4. Runs GRN inference for both canonical (gene-based) and AS-aware (transcript-based) inputs.
    5. Performs isoform categorization for downstream filtering.
    6. Aggregates and filters the inferred networks based on importance and frequency thresholds.
    7. Maps edges for comparison and categorizes them as gene-exclusive, isoform-exclusive, or common.
    8. Filters isoform-exclusive edges based on transcript plausibility.

    Parameters:
        config (dict): Configuration dictionary with file paths and parameters.
        transcript_tfs (pd.DataFrame): Transcript-level expression data for transcription factors.
            Columns must be sample IDs, rows must be Ensembl transcript IDs.
        gene_tfs (pd.DataFrame): Gene-level expression data for transcription factors.
            Columns must be sample IDs, rows must be Ensembl gene IDs.
        targets (pd.DataFrame): Expression data of target genes.
            Columns must be sample IDs, rows must be Ensembl gene IDs.


    Returns:
        pd.DataFrame: Plausibility-filtered edges found only in the isoform-aware regulatory network.
 
    '''
    
    # important files to load
    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)
    
    # create TF Database
    tf_database = utils_database.create_transcipt_annotation_database(tf_list=tf_list, appris_path= config['appris'], digger_path=config['digger'])


    print('Prepare for inference')

    data_canonical, data_asware, target_gene_list = prepare_for_inference(transcript_tfs, gene_tfs, targets)

    # get isoform gene mapping
    
    isoform_categories = utils_network.isoform_categorization(transcript_tfs, gene_tfs)

    # Inference of GRN
    as_aware_grn, canonical_grn = inference(config, config['nruns'], tf_list=tf_list, data_canonical=data_canonical, data_asware=data_asware, target_gene_list=target_gene_list, aggregate=True)

    #filter aggregate
    as_aware_grn = utils_network.filter_aggregate(as_aware_grn, threshold_frequency=1, threshold_importance=0.3)
    canonical_grn = utils_network.filter_aggregate(canonical_grn, threshold_frequency=1, threshold_importance=0.3)

    net_AS = utils_network.add_edge_key(as_aware_grn, biomart, source_column = 'source')
    net_canonical = utils_network.add_edge_key(canonical_grn, biomart, type='canonical', source_column='source')

    # Categorize edges into gene-unique, isoform-unique, common
    common_edges = utils_network.get_common_edges(net_canonical, net_AS)
    gene_unique, isoform_unique = utils_network.get_diff(net_canonical, net_AS)

    print('Number of edges in each category')
    print('Number of edges in gene-exclusive: ', len(gene_unique))
    print('Number of edges in isoform-exclusive: ', len(isoform_unique))
    print('Number of edges in common interactions: ', len(common_edges))
    
    # plausibility filtering
    plausibility_filtered = utils_network.plausibility_filtering(config, isoform_unique, isoform_categories, tf_database)

    return plausibility_filtered



def main():
    parser = argparse.ArgumentParser(description="Process a file from the command line.")
    
    # Add the file argument
    parser.add_argument('-f', type=str, help='Config File')
  
    # load config file
    args = parser.parse_args()
    with open(args.f, 'r') as f:
        config = yaml.safe_load(f)

    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)

    print('Load data')
    transcript_tfs, gene_tfs, targets = load_data(config, biomart, tf_list)
    
    inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets)
    
if __name__ == "__main__":
    main()
    