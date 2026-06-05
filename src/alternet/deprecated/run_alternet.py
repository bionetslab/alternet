import sys
sys.path.append('/data/bionets/og86asub/alternet-project/alternet/')
import pandas as pd
import alternet.data_preprocessing as preprocessing
from alternet.annotation import map_tf_ids
import alternet.annotation as annotation
import alternet.postprocessing as postprocessing
from alternet.inference import *
import yaml
import time
from typing import Any, Dict
import numpy as np

import seaborn as sns
import os.path as op
import os

def write_dict_to_yaml(data: Dict[str, Any], filepath: str):
    """Writes a Python dictionary to a YAML file."""

    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def add_gene_names_and_save(edge_data, prefix, net_type, filepath, biomart):
    outfile = op.join(filepath, f"{prefix}_{net_type}.tsv" )
    os.makedirs(filepath, exist_ok = True)

    # rename columns to fit a schema
    if 'target' in edge_data.columns and 'target_gene' not in edge_data.columns:
        edge_data = edge_data.rename(columns = {'target': 'target_gene'})


    if 'source_transcript' in edge_data.columns and 'source_gene_name' not in edge_data.columns:
        edge_data = postprocessing.add_transcript_names(edge_data, biomart, transcript_column = 'source_transcript')
    if 'source_gene' in edge_data.columns and 'source_gene_name' not in edge_data.columns:
        edge_data = postprocessing.add_gene_names(edge_data, biomart, gene_column = 'source_gene')
    if 'target_gene' in edge_data.columns and 'target_gene_name' not in edge_data.columns:
        edge_data = postprocessing.add_gene_names(edge_data, biomart, gene_column = 'target_gene')
        
    edge_data.to_csv(outfile)





def alternet_pipeline(transcript_data_cp, appris_df, digger_df, tf_list, biomart, results_path, prefix = 'condition', runs=10):

    tf_list = map_tf_ids(tf_list, biomart)
    tf_database = annotation.create_transcipt_annotation_database(tf_list=tf_list, appris_df= appris_df , digger=digger_df)
    transcript_mapper = annotation.create_transcript_mapping(biomart)

    
    gene_data_cp = transcript_data_cp.groupby('gene_id').sum()
    gene_data_cp = gene_data_cp.drop(columns={'transcript_id'})

    gene_data_cp.index.name = 'gene_id'
    gene_data_cp = gene_data_cp.reset_index()

    gene_target_names = list(tf_list['Gene stable ID'].unique())


    gene_data_cp = gene_data_cp.set_index('gene_id')
    gene_data_cp = gene_data_cp.T

    # scale data!!
    gene_data_cp_scaled = standardize_dataframe(gene_data_cp)

    transcript_data_cp = transcript_data_cp.set_index('transcript_id')
    transcript_data_cp = transcript_data_cp.drop('gene_id', axis=1)
    transcript_data_cp = transcript_data_cp.T

    # scale data!!
    transcript_data_cp_scaled = standardize_dataframe(transcript_data_cp)

    # compute isoform categories
    isoform_categories = postprocessing.isoform_categorization(transcript_data_cp, gene_data_cp, tf_list)
    gene_categories = postprocessing.get_gene_cases(isoform_categories)
    gene_to_transcript_mapping = annotation.create_filtered_gene_to_transcripts_mapping(biomart, gene_list = gene_data_cp_scaled.columns, transcript_list = transcript_data_cp_scaled.columns)


    start = time.monotonic()
    canonical_grn = inference(gene_data=gene_data_cp, tf_list=gene_target_names, target_names = 'all', n_runs=runs)
    end = time.monotonic()
    runtime = {'canonical': end-start}
    
    canonical_grn.to_csv(op.join(results_path, f"{prefix}_canonical.tsv"), header = True)

    hybrid_data = create_hybrid_data(transcript_data_cp, gene_data_cp, tf_list)
    hybrid_tf_names = list(tf_list['Transcript stable ID'].unique())
    target_names = list(gene_data_cp.columns)

    start = time.monotonic()
    as_aware_grn = inference(gene_data=hybrid_data, tf_list=hybrid_tf_names, target_names=target_names, n_runs = runs)
    runtime['as_aware'] = time.monotonic()-start

    as_aware_grn.to_csv(op.join(results_path, f"{prefix}_as_aware.tsv"), header = True) 
    write_dict_to_yaml(runtime, op.join(results_path, f"{prefix}_runtime.yaml"))

    start_postprocessing = time.monotonic()


    as_aware_grn = postprocessing.map_transcript_to_gene(as_aware_grn, transcript_mapper)
    as_aware_grn = postprocessing.create_edge_key(as_aware_grn)
    canonical_grn = postprocessing.create_edge_key(canonical_grn, source_column='source', target_column='target')
    canonical_grn = canonical_grn.rename(columns={'source': 'source_gene'})


    # Categorize edges into gene-unique, isoform-unique, common
    common_edges = postprocessing.get_common_edges(canonical_grn, as_aware_grn)
    gene_unique, isoform_unique = postprocessing.get_diff(canonical_grn, as_aware_grn)

    filtering_tracker = {'as_aware_base': as_aware_grn.shape[0], 'canonical_base': canonical_grn.shape[0],   'gene_unique_base': gene_unique.shape[0], 'isoform_unique_base': isoform_unique.shape[0]}


    # remove low quality edges
    gene_unique, gene_unique_filter_1 = postprocessing.frequency_filter(gene_unique,  threshold_frequency=runs)
    isoform_unique, isoform_unique_filter_1 = postprocessing.frequency_filter(isoform_unique, threshold_frequency=runs)

    # remove implausible edges, i.e. isoform edges where there is a dominant/single isoform which
    # should have been in gene GRN as well.
    isoform_unique, isoform_unique_plausibility_filter  = postprocessing.plausibility_filtering(isoform_unique, isoform_categories)
    # Merge isoform exon/domain information and save
    isoform_unique = annotation.annotate_isoform_exclusive_edges(isoform_unique, tf_database, transcript_column='source_transcript')
    


    # remove implausible edges, i.e. edges where there is a single/dominant isoform which should have been in
    # isoform GRN as well.
    gene_unique, gene_unique_plausibility_filter = postprocessing.plausibility_filtering_gene_unique(gene_unique, gene_categories)
    gene_unique = annotation.annotate_gene_exclusive_edges(gene_unique, annotation_database=tf_database, gene_transcript_mapping=gene_to_transcript_mapping)


    # remainder is in both GRNs. Merge the gene information to the isoform information
    merged_edges = postprocessing.create_common_edge_dataframe(common_edges)
    filtering_tracker['common_base'] = merged_edges.shape[0]

    # Split into dominant/single isoform and ambigous (multiple isoforms)
    consistent, ambigous = postprocessing.split_by_isoform_category(merged_edges, gene_categories)
    # remove edges which have vastly different frequencies and importances from dataframe (edges in common dataframe where there is only one single or
    # a dominant isoform should have a consistent importance and frequency).
    consistent, consistent_freq_filer = postprocessing.frequency_filtering_common_edges_dominant(consistent, threshold_frequency = runs)
    consistent, common_dominant_filter = postprocessing.plausibility_filtering_common_edges_dominant(consistent)
    consistent = annotation.annotate_consistent_edges(consistent, tf_database, transcript_column='source_transcript')


    # From remaining ambigous dataframe, select those edges which are sig. more importantn in isoform data frame
    likely_isoform_specific, ambigous, common_isoform_likely_filter = postprocessing.find_likely_isoform_specific(ambigous)
    # add transcript annotation and save
    likely_isoform_specific = annotation.annotate_isoform_exclusive_edges(likely_isoform_specific, tf_database, transcript_column='source_transcript')

    
    # from remaining ambigous dataframe, select those where the gene edge is more important than any isoform edge.
    likely_gene_specific, ambigous, common_gene_likely_filter = postprocessing.find_likely_gene_specific(ambigous)
    likely_gene_specific = annotation.annotate_gene_exclusive_edges(likely_gene_specific, annotation_database=tf_database, gene_transcript_mapping=gene_to_transcript_mapping)
    
    # remaining edges are ambigous.
    ambigous, ambigous_freq_filter = postprocessing.frequency_filtering_common_edges_dominant(ambigous, threshold_frequency = runs)

    correlation_collector = annotation.compute_isoform_gene_correlations(transcript_data_cp_scaled, gene_data_cp_scaled, gene_to_transcript_mapping)

    # Compute absolute importance threshold based on all surviving edgeds
    all_importances = np.concatenate([isoform_unique.mean_importance, gene_unique.median_importance, likely_gene_specific.median_importance_ge, likely_isoform_specific.median_importance_te, consistent.median_importance_te, ambigous.median_importance_te])
    absolute_threshold = np.percentile(all_importances, q=80)

    isoform_unique, filter_info_iu = postprocessing.filter_importance(isoform_unique, absolute_treshold=absolute_threshold,importance_column='median_importance' )
    gene_unique, filter_info_gu = postprocessing.filter_importance(gene_unique, absolute_treshold=absolute_threshold,importance_column='median_importance' )
    likely_gene_specific, filter_info_lgu =  postprocessing.filter_importance(likely_gene_specific, absolute_treshold=absolute_threshold,importance_column='median_importance_ge' )
    likely_isoform_specific, filter_info_liu =  postprocessing.filter_importance(likely_isoform_specific, absolute_treshold=absolute_threshold,importance_column='median_importance_ge' )
    consistent, filter_info_c  =  postprocessing.filter_importance(consistent, absolute_treshold=absolute_threshold,importance_column='median_importance_te' )
    ambigous, filter_info_a  =  postprocessing.filter_importance(ambigous, absolute_treshold=absolute_threshold,importance_column='median_importance_te' )


    add_gene_names_and_save(isoform_unique, prefix, 'unique_isoforms', results_path, biomart)
    add_gene_names_and_save(gene_unique, prefix, 'unique_genes', results_path,biomart)
    add_gene_names_and_save(consistent, prefix, 'consistent_both', results_path,biomart)
    add_gene_names_and_save(likely_isoform_specific, prefix, 'likely_isoform_specific', results_path,biomart)
    add_gene_names_and_save(likely_gene_specific, prefix, 'likely_gene_specific', results_path,biomart)
    add_gene_names_and_save(ambigous, prefix, 'ambigous', results_path,biomart)



    filtering_tracker['gene_unique_filter_frequency'] = gene_unique_filter_1
    filtering_tracker['isoform_unique_filter_frequency'] = isoform_unique_filter_1
    filtering_tracker['isoform_unique_plausibility_filter'] = isoform_unique_plausibility_filter
    filtering_tracker['gene_unique_plausibility_filter'] = gene_unique_plausibility_filter

    filtering_tracker['common_consistent_frequency_filter'] = consistent_freq_filer
    filtering_tracker['common_consitent_dominant_filter'] = common_dominant_filter
    
    filtering_tracker['common_isoform_likely_filter'] = common_isoform_likely_filter
    filtering_tracker['common_gene_likely_filter'] = common_gene_likely_filter
    filtering_tracker['common_ambigous_frequency_filter'] = ambigous_freq_filter

    filtering_tracker['isoform_unique_importance_threshold'] = filter_info_iu
    filtering_tracker['gene_unique_importance_threshold'] = filter_info_gu
    filtering_tracker['likely_gene_specific_importance_threshold'] = filter_info_lgu
    filtering_tracker['likely_isofrom_specific_importance_threshold'] = filter_info_liu
    filtering_tracker['common_consistent_importance_threshold'] = filter_info_c
    filtering_tracker['common_ambigous_importance_threshold'] = filter_info_a
    
    write_dict_to_yaml(filtering_tracker, op.join(results_path, f"{prefix}_filtering_tracker.yaml"))
    elapsed = {'postprocessing_elapsed': time.monotonic() - start_postprocessing}
    write_dict_to_yaml(elapsed, op.join(results_path, f"{prefix}_runtime_postprocessing.yaml"))





if __name__ == '__main__':

    results_path = "/data/bionets/og86asub/alternet-project/alternet-manuscript/results_new/"
    data_path = "/data/bionets/og86asub/alternet-project/alternet/data/"
    appris_path = "/data/bionets/og86asub/alternet-project/alternet/data/appris_data.appris.txt"
    digger_path = "/data/bionets/og86asub/alternet-project/alternet/data/digger_data.csv"
    biomart_path = '/data/bionets/og86asub/alternet-project/alternet/data/biomart.txt'
    tf_list_path = '/data/bionets/og86asub/alternet-project/alternet/data/allTFs_hg38.txt'


    os.makedirs(results_path, exist_ok=True)

    appris_df = pd.read_csv(appris_path, sep='\t')
    biomart = pd.read_csv(biomart_path, sep='\t')
    tf_list = pd.read_csv(tf_list_path, sep='\t', header = None)
   
    magnet_path = op.join(data_path, f"{prefix}_magnet_prefiltered_tpm.tsv")
    transcript_data_cp = pd.read_csv(magnet_path, sep='\t')
    