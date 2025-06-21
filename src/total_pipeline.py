import pandas as pd
from inference import *
from load_data import *
import utils_database
import utils_network
import json
import yaml
import argparse



def inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets, nruns):

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
    parser.add_argument('-n', type=int, help='Number of runs of grnboost')

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
    