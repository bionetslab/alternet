
import time
from alternet.inference import inference

import yaml
def write_dict_to_yaml(data, filepath):
    """Write dictionary to YAML file."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

        

runtime = {}
start = time.monotonic()
canonical_grn = inference(
    gene_data=gene_data_scaled,
    tf_list=tf_genes_in_data,
    target_names='all',
    n_runs=N_RUNS
)
runtime['canonical'] = time.monotonic() - start


canonical_grn = canonical_grn.rename(columns={'source': 'source_gene', 'target': 'target_gene'})
canonical_grn['reg_type'] = canonical_grn['source_gene'].map(gene_to_regtype)

canonical_grn.to_csv(op.join(results_path, f"{CONDITION}_canonical_raw.tsv"), sep='\t', index=False)






# Create hybrid data (TF transcripts + target genes)
hybrid_data = create_hybrid_data(
    transcript_data_matrix,  
    gene_data_matrix,        
    tf_list
)

start = time.monotonic()
as_source_grn = inference(
    gene_data=hybrid_data,
    tf_list=tf_transcripts_in_data,
    target_names=target_genes,
    n_runs=N_RUNS
)
runtime['as_aware_source'] = time.monotonic() - start


as_source_grn = as_source_grn.rename(columns={'source': 'source_transcript', 'target': 'target_gene'})
as_source_grn['source_gene'] = as_source_grn['source_transcript'].map(tx2gene)
as_source_grn['reg_type'] = as_source_grn['source_transcript'].map(tx_to_regtype)

as_source_grn.to_csv(op.join(results_path, f"{CONDITION}_as_aware_source_raw.tsv"), sep='\t', index=False)



start = time.monotonic()
fully_as_grn = inference(
    gene_data=transcript_data_scaled,
    tf_list=regulator_transcripts_in_data,
    target_names='all',
    n_runs=N_RUNS
)
runtime['fully_as_aware'] = time.monotonic() - start


# rename columns
fully_as_grn = fully_as_grn.rename(columns={'source': 'source_transcript', 'target': 'target_transcript'})
fully_as_grn['source_gene'] = fully_as_grn['source_transcript'].map(tx2gene)
fully_as_grn['target_gene'] = fully_as_grn['target_transcript'].map(tx2gene)
fully_as_grn['reg_type'] = fully_as_grn['source_transcript'].map(tx_to_regtype)
fully_as_grn.to_csv(op.join(results_path, f"{CONDITION}_fully_as_aware_raw.tsv"), sep='\t', index=False)




runtime['total'] = runtime['canonical'] + runtime['as_aware_source'] + runtime['fully_as_aware']
write_dict_to_yaml(runtime, op.join(results_path, f"{CONDITION}_runtime.yaml"))
print(f"\nTotal inference time: {runtime['total']/60:.2f} minutes")
