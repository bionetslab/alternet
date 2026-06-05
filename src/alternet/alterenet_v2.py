import numpy as np
import os
import os.path as op
import pandas as pd
from alternet.splicefactor_evidence import *
from alternet.compare_nets import *

from alternet.gtex_dataloader import *
from alternet.annotation import *




gtex_data_dir = '/data/bionets/datasets/hackathon/data/GTEX'
tissues = ['Blood', 'Brain', 'Adipose Tissue', 'Muscle', 'Blood Vessel','Heart']
# tissues = [ 'Ovary', 'Uterus', 'Vagina', 'Breast', 'Skin',
#        'Salivary Gland', 'Adrenal Gland', 'Thyroid', 'Lung', 'Spleen',
#        'Pancreas', 'Esophagus', 'Stomach', 'Colon', 'Small Intestine',
#        'Prostate', 'Testis', 'Nerve', 'Pituitary', 'Liver', 'Kidney',
#        'Cervix Uteri', 'Fallopian Tube', 'Bladder', 'Bone Marrow']


def alternet_main():

    data_path = "/data/bionets/og86asub/alternet-project/alternet/data"
    results_path = "/data/bionets/og86asub/alternet-project/alternet/raw_networks/"

    results_path_newa = "/data/bionets/og86asub/alternet-project/alternet/results_v2/"
    os.makedirs(results_path_newa, exist_ok = True)


    # Reference files
    appris_path = "appris_data.appris.txt"
    digger_path = "digger_data.csv"
    biomart_path = "biomart.txt"
    tf_list_path = "allTFs_hg38.txt"
    sf_list_path = "splicefactors.csv"

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



    for TISSUE in tissues:
        print(TISSUE)
        try:
            if op.exists(op.join(results_path_newa, f'{TISSUE}.tf.tsv')):
                print(f'File exisits {TISSUE}')
                continue
            
            transcript_data = pd.read_csv(op.join(gtex_data_dir, f'{TISSUE}.tsv'), sep = '\t', index_col = 0)

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
            canonical_grn = canonical_grn.rename(columns = {'source_gene': 'source', 'target_gene': 'target'})
            canonical_grn = canonical_names(canonical_grn, tx2gene, gene2tx)
            canonical_grn = filter_edges(canonical_grn)



            networks = pd.concat([canonical_grn, as_source_grn, fully_as_aware])
            networks = get_best_variable(networks)



            gene_dominance, gene_n_isoforms, tx_expression_share  = compute_dominance_metrics(transcript_data, sample_cols)
            networks = plausibility_filtering(networks, gene_dominance, r_dom = 0.9)

            usage_df, reliability_df = calculate_transcript_usage(transcript_data)
            transcript_data_temp = transcript_data.set_index('transcript_id')[sample_cols]
            usage_df = usage_df.set_index('transcript_id')[sample_cols]

            tf_net = networks[(networks.reg_type == 'TF')].copy()
            # Decide if TF_SF is splice factor or transcription factor
            nn = tf_sf_disambigouation_fully_as_aware(networks[(networks.reg_type == 'TF_SF') & (networks.target_type=='transcript')], regulator_list, transcript_data_temp, usage_df, reliability_df,
                                                    RHO_MIN = 0.3,  Q_MIN = 0.05, DU_MIN = 0.1, n_cores = 16)


            ab = compute_set_c(nn[nn.tfsf_category == 'tfsf_sf_like'], transcript_data, gene2tx, usage_df, reliability_df, sample_cols, epsilon=1e-6, n_cores=16)
            ac = compute_set_c(networks[(networks.reg_type == 'SF') & (networks.target_type=='transcript') & (networks.source_type=='transcript')], transcript_data, gene2tx, usage_df, reliability_df, sample_cols, epsilon=1e-6, n_cores=16)
            sf_net = pd.concat([ab, ac])


            # Concatenate TF like regulators
            ad = nn[nn.tfsf_category == 'tfsf_tf_like']
            tf_net = pd.concat([ad, tf_net])
            ambi_net = nn[nn.tfsf_category.isin(['tfsf_joint', 'tfsf_ambiguous'])]

            sf_net.to_csv(op.join(results_path_newa, f'{TISSUE}.sf.tsv'), sep='\t')
            ambi_net.to_csv(op.join(results_path_newa, f'{TISSUE}.ambi.tsv'), sep='\t')
            tf_net.to_csv(op.join(results_path_newa, f'{TISSUE}.tf.tsv'), sep='\t')

        except:
            print('Failure!')
    

if __name__ == '__main__':
    alternet_main()