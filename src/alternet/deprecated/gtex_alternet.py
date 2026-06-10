import numpy
from alternet.annotation import *
from alternet.data_preprocessing import standardize_dataframe
import numpy
import os
import os.path as op

from alternet.gtex_dataloader import *

from alternet import postprocessing
from alternet.deprecated.alternet_class import *
from alternet.edge_categorization import *



def write_dict_to_yaml(data, filepath):
    """Write dictionary to YAML file."""
    with open(filepath, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)


def main():
    
    

    data_path = "/data/bionets/og86asub/alternet-project/alternet/data"
    network_path = "/data/bionets/og86asub/alternet-project/alternet/raw_networks/"
    results_path = "/data/bionets/og86asub/alternet-project/alternet/results_gtex"

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
    TISSUE = "Heart"

    tissues = ['Blood', 'Brain', 'Adipose Tissue', 'Muscle', 'Blood Vessel',
       'Heart', 'Ovary', 'Uterus', 'Vagina', 'Breast', 'Skin',
       'Salivary Gland', 'Adrenal Gland', 'Thyroid', 'Lung', 'Spleen',
       'Pancreas', 'Esophagus', 'Stomach', 'Colon', 'Small Intestine',
       'Prostate', 'Testis', 'Nerve', 'Pituitary', 'Liver', 'Kidney',
       'Cervix Uteri', 'Fallopian Tube', 'Bladder', 'Bone Marrow']

    for TISSUE in tissues:

        try:
            CONDITION = TISSUE
            print(CONDITION)

            os.makedirs(results_path, exist_ok=True)

            biomart = pd.read_csv(op.join(data_path, biomart_path), sep='\t')
            tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
            gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()
            appris_df = pd.read_csv(op.join(data_path,appris_path), sep='\t')
            digger_df = pd.read_csv(op.join(data_path,digger_path), low_memory=False)
            # Load and map TF list
            tf_list_raw = pd.read_csv(op.join(data_path,tf_list_path), sep='\t', header=None)
            tf_list = map_tf_ids(tf_list_raw, biomart)


            # Load and map SF list
            sf_list_raw = pd.read_csv(op.join(data_path, sf_list_path), header=0, sep = ',')
            sf_list = map_sf_ids(sf_list_raw.loc[:, ['Splicing_Factor']], biomart)
            # Combine TF and SF lists
            regulator_list = combine_tf_sf_lists(tf_list, sf_list)


            tx_to_regtype = dict(zip(regulator_list['Transcript stable ID'], regulator_list['Regulator_type']))
            gene_to_regtype = regulator_list.groupby('Gene stable ID')['Regulator_type'].first().to_dict()
            gtex_data_dir = '/data/bionets/datasets/hackathon/data/GTEX'
            params = {'sample_attributes': op.join(gtex_data_dir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), 'tissue': CONDITION, 'transcript_data':op.join(gtex_data_dir, 'GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct')}
            tissue_ids = retrieve_GTEX_tissue_sampleids(params['sample_attributes'], tissue=params['tissue'])
            transcript_data = read_GTEX_transcript_expression(params['transcript_data'], tissue_ids)
            transcript_data = clean_GTEX_tissue_transcript_counts(transcript_data, biomart)
            transcript_data = variance_filtering(transcript_data)


            canonical_grn = pd.read_csv(op.join(network_path, f"{CONDITION}_canonical_raw.tsv"), sep='\t')
            as_source_grn = pd.read_csv(op.join(network_path,  f"{CONDITION}_as_aware_source_raw.tsv"), sep='\t')
            fully_as_grn= pd.read_csv(op.join(network_path, f"{CONDITION}_fully_as_aware_raw.tsv"), sep='\t')


            alternet_obj21 = Alternet(canonical_grn, as_source_grn, fully_as_grn, transcript_data, regulator_list, tx_to_regtype , 'gene_id', 'transcript_id', min_frequency=10, importance_percentile = 0.8)


            os.makedirs(op.join(results_path, CONDITION), exist_ok = True)
            alternet_obj21.set_a.to_csv(op.join(results_path, CONDITION, f"{CONDITION}_set_a.tsv"), sep = '\t')
            alternet_obj21.set_b.to_csv(op.join(results_path, CONDITION, f"{CONDITION}_set_b.tsv"), sep = '\t')
            alternet_obj21.set_c.to_csv(op.join(results_path, CONDITION, f"{CONDITION}_set_c.tsv"), sep = '\t')
            alternet_obj21.set_d.to_csv(op.join(results_path, CONDITION, f"{CONDITION}_set_d.tsv"), sep = '\t')
        
        except:
            print(f'Failure on tissue: {CONDITION}')


if __name__ == "__main__":
    main()