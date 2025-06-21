import pandas as pd
from inference import *
from load_data import *
import utils_database
import utils_network
import json
import yaml
import argparse


def main():
    parser = argparse.ArgumentParser(description="Process a file from the command line.")
    
    # Add the file argument
    parser.add_argument('-f', type=str, help='Config File')
    parser.add_argument('-n', type=int, help='Number of runs of grnboost')

    # load config file
    args = parser.parse_args()
    with open(args.f, 'r') as f:
        config = yaml.safe_load(f)

    # important paths, prefixes, files:
    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')

    as_aware_prefix = f"{config['tissue']}_as-aware.network_"
    canonical_prefix = f"{config['tissue']}_canonical.network_"

    as_aggregated = op.join(results_dir_grn, as_aware_prefix + 'aggregated.tsv')
    canonical_aggregated = op.join(results_dir_grn, canonical_prefix + 'aggregated.tsv')


    # important files to load
    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)
    

    # create TF Database
    tf_database = utils_database.create_transcipt_annotation_database(tf_list=tf_list,appris_path= config['appris'], biomart_digger_mapping=config['biomart_digger_mapping'], digger_path=config['digger'])


    print('Load data and prepare for inference')
    
    transcript_tfs, gene_tfs, targets = load_data(config, biomart, tf_list)
    data_canonical, data_asware, target_gene_list = prepare_for_inference(transcript_tfs, gene_tfs, targets)

    # get isoform gene mapping
    
    isoform_categories = utils_network.isoform_categorization(transcript_tfs, gene_tfs)


    # Inference of GRN
    as_aware_grn, canonical_grn = inference(config, args.n, tf_list=tf_list, data_canonical=data_canonical, data_asware=data_asware, target_gene_list=target_gene_list, aggregate=True)

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
    
    AS_only = isoform_unique.merge(isoform_categories.loc[:, ['transcript_id', 'isoform_category']], left_on='source', right_on='transcript_id', how='left').drop(columns=['transcript_id'])

    AS_only_annotated = utils_database.merge_annotations_to_grn(AS_only, tf_database)

    file = op.join(results_dir_grn, as_aware_prefix + f"annotated_w_database.tsv")
    AS_only_annotated.to_csv(file, sep='\t', index=False)

    #filter
    protein_coding = AS_only_annotated[AS_only_annotated['Protein Coding'] == True]
    protein_coding = protein_coding.drop(columns=['mean_importance', 'source_transcript', 'target_gene', 'type', 'edge_key'])
    protein_coding = protein_coding[protein_coding['isoform_category'] != 'dominant']
    sort_pc = protein_coding.sort_values(by=['frequency', 'median_importance'], ascending=[False, False])
    
    file = op.join(results_dir_grn, as_aware_prefix + f"plausibility_filtered_iso_unique.tsv")
    sort_pc.to_csv(file, sep='\t', index=False)


if __name__ == "__main__":
    main()
    