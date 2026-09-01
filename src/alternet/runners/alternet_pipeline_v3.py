import numpy as np
import os
import os.path as op
import pandas as pd
import time
import yaml
import argparse


from alternet.splicefactor_evidence import compute_dominance_metrics, calculate_transcript_usage, tf_sf_disambigouation_fully_as_aware, compute_set_c
from alternet.compare_nets import *
from alternet.annotation import *
from alternet.data_preprocessing import *



def write_dict_to_yaml(data, filepath):
    """Write dictionary to YAML file."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def main():
    parser = argparse.ArgumentParser(
        description="Run alternet downstream"
    )

    parser.add_argument(
        "--transcript_file", 
        type=str, 
        required=True, 
        help="Path to the directory containing GTEx raw data."
    )

    parser.add_argument(
        "--canonical", 
        type=str, 
        required=True, 
        help="Path to the directory containing GTEx raw data."
    )
    
    parser.add_argument(
        "--source_as", 
        type=str, 
        required=True, 
        help="Path to the directory containing GTEx raw data."
    )
    
    parser.add_argument(
        "--fully_as", 
        type=str, 
        required=True, 
        help="Path to the directory containing GTEx raw data."
    )

    parser.add_argument(
        "--results_path", 
        type=str, 
        required=True, 
        help="Path where the output results will be saved."
    )

    parser.add_argument(
        "--experiment_name", 
        type=str, 
        required=True, 
        help="The specific tissue type to analyze (e.g., 'Liver', 'Brain')."
    )

    parser.add_argument(
        "--p_value", 
        type=float, 
        required=False, 
        default=1,
        help="Minimal required number of edge occurences"
    )

    parser.add_argument(
        "--biomart_file", 
        type=str, 
        help="Minimal required number of edge occurences"
    )


    parser.add_argument(
        "--regulator_file", 
        type=str, 
        help="Minimal required number of edge occurences"
    )

    parser.add_argument(
        "--appris", 
        type=str, 
        help="Minimal required number of edge occurences"
    )

    parser.add_argument(
        "--digger", 
        type=str, 
        help="Minimal required number of edge occurences"
    )



    # 3. Parse the arguments
    args = parser.parse_args()

    results_path = args.results_path
    experiment_name = args.experiment_name
    
    p_value = args.p_value

    os.makedirs(results_path, exist_ok = True)

    regulator_list  = pd.read_csv(args.regulator_file, sep = '\t')
    biomart = pd.read_csv(args.biomart_file, sep='\t')
    tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
    gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()
    gene2regtype = dict(zip(regulator_list['Transcript stable ID'], regulator_list['Regulator_type'])) | dict(zip(regulator_list['Gene stable ID'], regulator_list['Regulator_type']))

   
    appris_df = pd.read_csv(args.appris, sep='\t')
    digger_df = pd.read_csv(args.digger, low_memory=False)

    tf_database = create_transcipt_annotation_database(
        tf_list=regulator_list, appris_df=appris_df, digger=digger_df
    )

    transcript_data = pd.read_csv(args.transcript_file, sep = '\t', index_col = 0)
    sample_cols = [c for c in transcript_data.columns if c not in ['transcript_id', 'gene_id']]



    as_source_grn = pd.read_csv(args.source_as, sep='\t')
    print(as_source_grn)
    as_source_grn = canonical_names(as_source_grn, tx2gene, gene2tx)
    as_source_grn = filter_edges(as_source_grn, frequency=10, imp_col = 'median_importance')

    fully_as_aware = pd.read_csv(args.fully_as, sep='\t')
    fully_as_aware = canonical_names(fully_as_aware, tx2gene, gene2tx)
    fully_as_aware = filter_edges(fully_as_aware,  frequency=10, imp_col = 'median_importance')


    canonical_grn = pd.read_csv(args.canonical, sep='\t')
    canonical_grn = canonical_names(canonical_grn, tx2gene, gene2tx)
    canonical_grn = filter_edges(canonical_grn, frequency=10, imp_col = 'median_importance')

    print(canonical_grn)

    networks = pd.concat([canonical_grn, as_source_grn, fully_as_aware])
    networks = get_best_variable(networks)



    gene_dominance, gene_n_isoforms, tx_expression_share  = compute_dominance_metrics(transcript_data, sample_cols)
    networks = plausibility_filtering(networks, gene_dominance, r_dom = 0.9)
    print(networks)

    usage_df, reliability_df = calculate_transcript_usage(transcript_data)
    transcript_data_temp = transcript_data.set_index('transcript_id')[sample_cols]
    usage_df = usage_df.set_index('transcript_id')[sample_cols]

    # Source == TF net
    tf_net = networks[(networks.reg_type == 'TF')].copy()

    # source == SF net
    sf_candidates = networks[(networks.reg_type == 'SF') & (networks.target_type=='transcript') & (networks.source_type=='transcript')]
    if sf_candidates.shape[0]>0:
        sf_net = compute_set_c(sf_candidates, transcript_data, gene2tx, usage_df, reliability_df, sample_cols, epsilon=1e-6, n_cores=16)
    else:
        sf_net = None

    print(sf_net)


    # Decide if TF_SF is splice factor or transcription factor
    if networks[(networks.reg_type == 'TF_SF')].shape[0]>0:
        nn = tf_sf_disambigouation_fully_as_aware(networks[(networks.reg_type == 'TF_SF') & (networks.target_type=='transcript')], regulator_list, transcript_data_temp, usage_df, reliability_df,
                                                RHO_MIN = 0.3,  Q_MIN = 0.05, DU_MIN = 0.1, n_cores = 16)

        ab = compute_set_c(nn[nn.tfsf_category == 'tfsf_sf_like'], transcript_data, gene2tx, usage_df, reliability_df, sample_cols, epsilon=1e-6, n_cores=16)
        if sf_net is not None and ab.shape[0]>0:
            sf_net = pd.concat([ab, sf_net])


        # Concatenate TF like regulators if there are any
        ad = nn[nn.tfsf_category == 'tfsf_tf_like']
        tf_net = pd.concat([ad, tf_net])

        # the rest
        ambi_net = nn[nn.tfsf_category.isin(['tfsf_joint', 'tfsf_ambiguous'])]
    else:
        ambi_net = None



    if tf_net.shape[0]>0:
        tf_net = annotate_isoform_exclusive_edges(tf_net, tf_database, transcript_column='source_transcript')
        tf_net = annotate_isoform_exclusive_edges(tf_net, tf_database, transcript_column='target_transcript', suffixes = ('_source', '_target'))
        tf_net.to_csv(op.join(results_path, f'{experiment_name}.tf.tsv'), sep='\t')

    if (ambi_net is not None) and (ambi_net.shape[0]>0):
        ambi_net = annotate_isoform_exclusive_edges(ambi_net, tf_database, transcript_column='source_transcript')
        ambi_net = annotate_isoform_exclusive_edges(ambi_net, tf_database, transcript_column='target_transcript', suffixes = ('_source', '_target'))
        ambi_net.to_csv(op.join(results_path, f'{experiment_name}.ambi.tsv'), sep='\t')

    if (sf_net is not None) and sf_net.shape[0]>0:
        sf_net = annotate_isoform_exclusive_edges(sf_net, tf_database, transcript_column='source_transcript')
        sf_net = annotate_isoform_exclusive_edges(sf_net, tf_database, transcript_column='target_transcript', suffixes = ('_source', '_target'))
        sf_net.to_csv(op.join(results_path, f'{experiment_name}.sf.tsv'), sep='\t')

    print(op.join(results_path, f'{experiment_name}.sf.tsv'))




if __name__ == "__main__":
    main()