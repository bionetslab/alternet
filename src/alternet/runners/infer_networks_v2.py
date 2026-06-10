
import time
from alternet.inference import inference
from alternet.data_preprocessing import *
from alternet.annotation import *
import os.path as op
import yaml
import os
from alternet.gtex_dataloader import *
import pandas as pd


import argparse
import os
import sys

def write_dict_to_yaml(data, filepath):
    """Write dictionary to YAML file."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def main():
    parser = argparse.ArgumentParser(
        description="Process GTEx data for a specific tissue type."
    )

    parser.add_argument(
        "--transcript_file", 
        type=str, 
        required=True, 
        help="Path to the directory containing GTEx raw data."
    )

    parser.add_argument(
        "--data_path", 
        type=str, 
        required=True, 
        help="Path to the input data file."
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
        "--n_runs", 
        type=int, 
        required=False, 
        default=10,
        help="Number of times to run GRNboost for"
    )
    
    parser.add_argument(
        "--biomart_file", 
        type=int, 
        required=False, 
        default=10,
        help="Minimal required number of edge occurences"
    )


    parser.add_argument(
        "--regulator_file", 
        type=int, 
        required=False, 
        default=10,
        help="Minimal required number of edge occurences"
    )

    parser.add_argument(
        "--canonical", 
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



    # 3. Parse the arguments
    args = parser.parse_args()

    results_path = args.results_path
    experiment_name = args.experiment_name
    N_RUNS = 10


    regulator_list  = pd.read_csv(args.regulator_file, sep = '\t')
    biomart = pd.read_csv(args.biomart_file, sep='\t')
    tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))

    tx_to_regtype = dict(zip(regulator_list['Transcript stable ID'], regulator_list['Regulator_type']))
    gene_to_regtype = regulator_list.groupby('Gene stable ID')['Regulator_type'].first().to_dict()
    


    runtime_filepath = op.join(results_path_tissue, f"{experiment_name}_runtime.yaml")

    if op.exists(runtime_filepath):
        with open(runtime_filepath, 'r') as f:
            runtime = yaml.safe_load(f) or {}  # 'or {}' handles empty files safely
    else:
        runtime = {}


    results_path_tissue = op.join(results_path, experiment_name)
    os.makedirs(results_path_tissue, exist_ok=True)
    
    transcript_data = pd.read_csv(args.transcript_file, sep = '\t', index_col = 0)
    sample_cols = [c for c in transcript_data.columns if c not in ['transcript_id', 'gene_id']]


    gene_data = transcript_data.groupby('gene_id')[sample_cols].sum().reset_index()

    
    gene_data_matrix = gene_data.set_index('gene_id')[sample_cols].T
    transcript_data_matrix = transcript_data.set_index('transcript_id')[sample_cols].T

    gene_data_scaled = standardize_dataframe(gene_data_matrix)
    transcript_data_scaled = standardize_dataframe(transcript_data_matrix)

    regulator_transcripts_all = list(set(regulator_list['Transcript stable ID']) & set(transcript_data_scaled.columns))
    regulator_genes_all = list(set(regulator_list['Gene stable ID']) & set(transcript_data_scaled.columns))


    runtime['n_transcripts'] = transcript_data.shape[0]
    runtime['n_genes'] = gene_data.shape[0]


    # GENE-all network
    # Hybrid data = Transcription factor isoform counts, target gene counts
    # Use target genes to get only iso-gene edges
    # use TF transcripts as regulators
    # Create hybrid data (TF transcripts + target genes)
    hybrid_data = pd.concat([transcript_data, gene_data], axis =1)

    target_genes = list(gene_data.columns)
    if not op.exists(args.as_source):
        print(f"File not found. Running GRN inference for {experiment_name}... using transcripts as regulators")
        start = time.monotonic()
        as_source_grn = inference(
            gene_data=hybrid_data,
            tf_list=regulator_genes_all,
            target_names='all',
            n_runs=N_RUNS
        )
        runtime['as_aware_source'] = time.monotonic() - start

        as_source_grn['source_type'] = 'transcript'
        as_source_grn['target_type'] = as_source_grn['target'].apply(lambda x: 'gene' if str(x).startswith('ENSG') else 'transcript')
        as_source_grn['source_gene'] = as_source_grn['source_transcript'].map(tx2gene)
        as_source_grn['reg_type'] = as_source_grn['source_transcript'].map(tx_to_regtype)
        as_source_grn.to_csv(args.as_source, sep='\t', index=False)
        write_dict_to_yaml(runtime, op.join(results_path_tissue, f"{experiment_name}_runtime.yaml"))


    #ISO-all-network
    # Use transcript data 
    # use tf and sf isoforms as regulators.
    if not op.exists(args.fully_as):
        print(f"File not found. Running GRN inference for {experiment_name}... using transcripts as regulators and targets")
        start = time.monotonic()
        as_source_grn = inference(
            gene_data=hybrid_data,
            tf_list=regulator_transcripts_all,
            target_names='all',
            n_runs=N_RUNS
        )
        runtime['as_aware_full'] = time.monotonic() - start

        as_source_grn['source_type'] = 'transcript'
        as_source_grn['target_type'] = as_source_grn['target'].apply(lambda x: 'gene' if str(x).startswith('ENSG') else 'transcript')
        as_source_grn['source_gene'] = as_source_grn['source_transcript'].map(tx2gene)
        as_source_grn['reg_type'] = as_source_grn['source_transcript'].map(tx_to_regtype)
        as_source_grn.to_csv(args.as_source, sep='\t', index=False)
        write_dict_to_yaml(runtime, op.join(results_path_tissue, f"{experiment_name}_runtime.yaml"))


if __name__ == "__main__":
    main()
