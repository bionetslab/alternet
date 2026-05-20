import numpy as np
import numpy as np
import pandas as pd
import os 
import os.path as op
from alternet.postprocessing import *


r_dom = 0.9
from alternet.plausibility_filtering import *

def plausibility_filtering(df,gene_dominance, r_dom = 0.9):
    df['reg_dominance'] = df['source_gene'].map(gene_dominance).fillna(0.5)
    df['target_dominance'] = df['target_gene'].map(gene_dominance).fillna(0.5)

        # Initialize filter columns
    df['is_plausible'] = True

    # If both isoforms make a large proportion of the total expression, we should
    # at least find three edges for all more or less identical interactions.
    mask_iso_specific = ~((df.reg_dominance>r_dom) & (df.target_dominance>r_dom) & df['is_unique_edge'])
    df['is_plausible'] = mask_iso_specific

    return df



def canonical_names(df, tx2gene, gene2tx):

    df["source_type"] = np.where(df["source"].str.startswith("ENSG"), "gene", "transcript")

    df["target_type"] = np.where(df["target"].str.startswith("ENSG"), "gene", "transcript")

    df["source_gene"] = np.where(df["source_type"] == 'gene', df["source"], df["source"].map(tx2gene))

    df["target_gene"] = np.where(df["target_type"] == 'gene', df["target"], df["target"].map(tx2gene))

    df["source_transcript"] = np.where(df["source_type"] == "transcript", df["source"], np.nan )
    
    
    df["target_transcript"] = np.where(df["target_type"] == "transcript", df["target"], np.nan)

    return df


def get_best_variable(df, r_iso=0.66,r_eq=0.8, r_iso_strict=0.5, eps=1e-6):
    df['edge_gg'] = df['source_gene'] + '_' + df['target_gene']
    df['edge_type'] = df['source_type'] + '-' + df['target_type']

    # 1. Find the maximum 'mean_importance' per group and map it to all elements
    df['max_mean_importance'] = df.groupby('edge_gg')['mean_importance'].transform('max')
    
    df["is_unique_edge"] = df.groupby("edge_gg")["edge_gg"].transform("count") == 1

    # 2. Compute the ratio of each element's importance to its group's maximum
    df['importance_ratio'] = df['mean_importance'] / df['max_mean_importance'] + eps

    is_max = df['mean_importance'] == df['max_mean_importance']

    df['max_other_ratio'] = (df['importance_ratio'].where(~is_max)
    .groupby(df['edge_gg'])
    .transform('max')
    .fillna(0)  # If a group only has 1 element, the max other ratio is 0
    )

        # 5. Define the conditions for each category
    # Specific: The max element is way ahead of the runner-up
    cond_specific_strict = is_max & (df['max_other_ratio'] < r_iso_strict)

    cond_specific = is_max & (df['max_other_ratio'] < r_iso)

    # Equivalent: There is no singular max because the runner-up is within the 1/r_eq window
    # This labels both the max and any other elements within that top window as 'equivalent'
    cond_equivalent = (df['importance_ratio'] >= r_eq) & (df['max_other_ratio'] >= r_eq)

    # 6. Apply categories using np.select
    conditions = [cond_specific_strict, cond_specific, cond_equivalent ]
    choices = ['specific_strict', 'specific', 'equivalent']

    df['category'] = np.select(conditions, choices, default='other')

    df['n_equivalent_edge_types'] = (df['edge_type']
        .where(cond_equivalent)
        .groupby(df['edge_gg'])
        .transform('nunique')
    )
    # Clean up helper columns
    df = df.drop(columns=['max_mean_importance', 'max_other_ratio'])
    return df


gtex_data_dir = '/data/bionets/datasets/hackathon/data/GTEX'

tissues = ['Blood', 'Brain', 'Adipose Tissue', 'Muscle', 'Blood Vessel',
       'Heart', 'Ovary', 'Uterus', 'Vagina', 'Breast', 'Skin',
       'Salivary Gland', 'Adrenal Gland', 'Thyroid', 'Lung', 'Spleen',
       'Pancreas', 'Esophagus', 'Stomach', 'Colon', 'Small Intestine',
       'Prostate', 'Testis', 'Nerve', 'Pituitary', 'Liver', 'Kidney',
       'Cervix Uteri', 'Fallopian Tube', 'Bladder', 'Bone Marrow']


def alternet_main():

    data_path = "/data/bionets/og86asub/alternet-project/alternet/data"
    results_path = "/data/bionets/og86asub/alternet-project/alternet/raw_networks/"

    results_path_newa = "/data/bionets/og86asub/alternet-project/alternet/results_v2/"
    os.makedirs(results_path_newa, exist_ok = True)


    # Reference files
    appris_path = "appris_data.appris.txt"
    digger_path = "digger_data.csv"
    biomart_path = "biomart.txt"


    biomart = pd.read_csv(op.join(data_path, biomart_path), sep='\t')
    tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
    gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()


    for TISSUE in tissues:
        try:
            transcript_data = pd.read_csv(op.join(gtex_data_dir, f'{TISSUE}.tsv'), sep = '\t')

            sample_cols = [c for c in transcript_data.columns if c not in ['transcript_id', 'gene_id']]

            as_source_grn = pd.read_csv(op.join(results_path,  f"{TISSUE}_as_aware_source_raw.tsv"), sep='\t')
            as_source_grn = as_source_grn.rename(columns = {'source_transcript':'source', 'target_gene': 'target'})
            as_source_grn = canonical_names(as_source_grn, tx2gene, gene2tx)
            as_source_grn = filter_edges(as_source_grn)
   
            fully_as_aware = pd.read_csv(op.join(results_path,  f"{TISSUE}_fully_as_aware_raw.tsv"), sep='\t')
            fully_as_aware = fully_as_aware.rename(columns = {'source_transcript':'source', 'target_transcript': 'target'})
            fully_as_aware = canonical_names(fully_as_aware, tx2gene, gene2tx)
            fully_as_aware = filter_edges(fully_as_aware)


            canonical_grn = pd.read_csv(op.join(results_path, f"{TISSUE}_canonical_raw.tsv"), sep='\t')
            canonical_grn = canonical_grn.rename(columns = {'source_gene': 'source'})
            canonical_grn = canonical_names(canonical_grn, tx2gene, gene2tx)
            canonical_grn = filter_edges(canonical_grn)


            networks = pd.concat([canonical_grn, as_source_grn, fully_as_aware])
            networks = get_best_variable(networks)

            gene_dominance, gene_n_isoforms, tx_expression_share  = compute_dominance_metrics(transcript_data, sample_cols)
            networks = plausibility_filtering(networks, gene_dominance, r_dom = 0.9)

            networks.to_csv(op.join(results_path_newa, f'{TISSUE}_alternetv3.tsv'), sep = '\t')
        except:
            print(f'Failure on tissue {TISSUE}')


if __name__ == '__main__':
    alternet_main()