
import time
from alternet.inference import inference
from alternet.data_preprocessing import *
import os.path as op
import yaml

if __name__ == '__main__':
    gtex_data_dir = '/data/bionets/datasets/hackathon/data/GTEX'
    

    data_path = "/data/bionets/og86asub/alternet-project/alternet/data"
    results_path = "/data/bionets/og86asub/alternet-project/alternet/results-2.0"

    # Reference files
    appris_path = "appris_data.appris.txt"
    digger_path = "digger_data.csv"
    biomart_path = "biomart.txt"
    tf_list_path = "allTFs_hg38.txt"
    sf_list_path = "splicefactors.csv"

    # Expression data
    gtex_transcript_tpm_path = "GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt"
    gtex_sample_attributes_path = "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"

    # Tissue to analyze
    TISSUE = "Bladder"

    # Number of GRNBoost2 runs
    N_RUNS = 10


    biomart = pd.read_csv(op.join(data_path, biomart_path), sep='\t')
    tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
    gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()

        # Load transcription factor list
    tf_list_raw = pd.read_csv(op.join(data_path, tf_list_path), sep='\t', header=None)
    tf_list = map_tf_ids(tf_list_raw, biomart)

    # Load and map SF list
    sf_list_raw = pd.read_csv(op.join(data_path, sf_list_path), header=0, sep = ',')
    sf_list = map_sf_ids(sf_list_raw.loc[:, ['Splicing_Factor']], biomart)

    # Combine TF and SF lists
    regulator_list = combine_tf_sf_lists(tf_list, sf_list)
    tx_to_regtype = dict(zip(regulator_list['Transcript stable ID'], regulator_list['Regulator_type']))
    gene_to_regtype = regulator_list.groupby('Gene stable ID')['Regulator_type'].first().to_dict()
    
    
    N_RUNs = 10

    for TISSUE in ['Bladder', 'Liver', 'Spleen']:
        
        CONDITION = TISSUE
        results_path_tissue = op.join(results_path, TISSUE)
        os.makedirs(results_path_tissue, exist_ok=True)
        
        params = {'sample_attributes': op.join(gtex_data_dir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), 'tissue': 'Liver', 'transcript_data':op.join(gtex_data_dir, 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct')}
        tissue_ids = retrieve_GTEX_tissue_sampleids(params['sample_attributes'], tissue=params['tissue'])
        transcript_data = read_GTEX_transcript_expression(params['transcript_data'], tissue_ids)
        transcript_data = clean_GTEX_tissue_transcript_counts(transcript_data, biomart)
        transcript_data = variance_filtering(transcript_data)


        tf_genes_in_data = list(set(tf_list['Gene stable ID']) & set(gene_data_scaled.columns))
        tf_transcripts_in_data = list(set(tf_list['Transcript stable ID']) & set(transcript_data_scaled.columns))
        regulator_transcripts_in_data = list(set(regulator_list['Transcript stable ID']) & set(transcript_data_scaled.columns))

        target_genes = list(gene_data_scaled.columns)
        target_transcripts = list(transcript_data_scaled.columns)


        sample_cols = [c for c in transcript_data.columns if c not in ['transcript_id', 'gene_id']]
        gene_data = transcript_data.groupby('gene_id')[sample_cols].sum().reset_index()

        # Create expression matrices (samples × features)
        gene_data_matrix = gene_data.set_index('gene_id')[sample_cols].T
        transcript_data_matrix = transcript_data.set_index('transcript_id')[sample_cols].T


        gene_data_scaled = standardize_dataframe(gene_data_matrix)
        transcript_data_scaled = standardize_dataframe(transcript_data_matrix)

        transcript_data_scaled, gene_data_scaled = remove_problematic_transcripts(transcript_data_scaled, gene_data_scaled)

        # Create hybrid data (TF transcripts + target genes)
        hybrid_data = create_hybrid_data(
            transcript_data_scaled,  
            gene_data_scaled,        
            tf_list
        )


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
        canonical_grn.to_csv(op.join(results_path_tissue, f"{CONDITION}_canonical_raw.tsv"), sep='\t', index=False)


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
        as_source_grn.to_csv(op.join(results_path_tissue, f"{CONDITION}_as_aware_source_raw.tsv"), sep='\t', index=False)




        start = time.monotonic()
        fully_as_grn = inference(
            gene_data=transcript_data_scaled,
            tf_list=regulator_transcripts_in_data,
            target_names='all',
            n_runs=N_RUNS
        )
        runtime['fully_as_aware'] = time.monotonic() - start

        fully_as_grn = fully_as_grn.rename(columns={'source': 'source_transcript', 'target': 'target_transcript'})
        fully_as_grn['source_gene'] = fully_as_grn['source_transcript'].map(tx2gene)
        fully_as_grn['target_gene'] = fully_as_grn['target_transcript'].map(tx2gene)
        fully_as_grn['reg_type'] = fully_as_grn['source_transcript'].map(tx_to_regtype)
        fully_as_grn.to_csv(op.join(results_path_tissue, f"{CONDITION}_fully_as_aware_raw.tsv"), sep='\t', index=False)



        def write_dict_to_yaml(data, filepath):
            """Write dictionary to YAML file."""
            with open(filepath, 'w') as f:
                yaml.dump(data, f, default_flow_style=False)

        runtime['total'] = runtime['canonical'] + runtime['as_aware_source'] + runtime['fully_as_aware']
        write_dict_to_yaml(runtime, op.join(results_path_tissue, f"{CONDITION}_runtime.yaml"))
        print(f"\nTotal inference time: {runtime['total']/60:.2f} minutes")
