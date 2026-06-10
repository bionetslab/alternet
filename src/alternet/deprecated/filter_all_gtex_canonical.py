import numpy
from alternet.annotation import *
from alternet.data_preprocessing import standardize_dataframe
import numpy
import os
import os.path as op

from alternet.gtex_dataloader import *

from alternet.postprocessing import filter_edges
import yaml


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

            canonical_grn = pd.read_csv(op.join(network_path, f"{CONDITION}_canonical_raw.tsv"), sep='\t')
            as_source_grn = pd.read_csv(op.join(network_path,  f"{CONDITION}_as_aware_source_raw.tsv"), sep='\t')
            fully_as_grn= pd.read_csv(op.join(network_path, f"{CONDITION}_fully_as_aware_raw.tsv"), sep='\t')

            canonical_grn = filter_edges(canonical_grn,  importance_percentile=0.9)
            as_source_grn = filter_edges(as_source_grn,  importance_percentile=0.9)
            fully_as_grn = filter_edges(fully_as_grn,  importance_percentile=0.9)

            canonical_grn.to_csv(op.join(network_path, f"{CONDITION}_canonical_filtered.tsv"), sep = '\t')
            as_source_grn.to_csv(op.join(network_path, f"{CONDITION}_as_aware_source_filtered.tsv"), sep = '\t')
            fully_as_grn.to_csv(op.join(network_path, f"{CONDITION}_fully_as_aware_filtered.tsv"), sep = '\t')


        except:
            print(f'Failure on tissue: {CONDITION}')


if __name__ == "__main__":
    main()