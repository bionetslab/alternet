import pandas as pd
from inference import *
from load_data import *
import os
import gseapy as gp
import magnet 


def load_network(config): 
    '''
    Loads aggregated isoform-level and gene-level GRNs.

    Parameters:
        config (dict): Configuration dictionary containing:
            - 'results_dir' (str): Base directory where results are stored.
            - 'tissue' (str): Tissue name used as a prefix to locate the GRN result files.

    Returns:
        tuple:
            - as_aware_grn (pd.DataFrame): Aggregated AS-aware GRN.
            - canonical_grn (pd.DataFrame): Aggregated canonical GRN.
    '''

    # important paths, prefixes, files:
    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')

    as_aware_prefix = f"{config['tissue']}_as-aware.network_"
    canonical_prefix = f"{config['tissue']}_canonical.network_"

    as_aggregated = op.join(results_dir_grn, as_aware_prefix + 'aggregated.tsv')
    canonical_aggregated = op.join(results_dir_grn, canonical_prefix + 'aggregated.tsv')

    as_aware_grn = pd.read_csv(as_aggregated, sep='\t')
    canonical_grn = pd.read_csv(canonical_aggregated, sep='\t')

    return as_aware_grn, canonical_grn

def load_plausibility_filtered(config):
    '''
    Loads the plausibility-filtered isoform-unique gene regulatory network (GRN).

    Parameters:
        config (dict): Configuration dictionary containing:
            - 'results_dir' (str): Base directory where results are stored.
            - 'tissue' (str): Tissue name used as a prefix to locate the GRN result files.

    Returns:
        pd.DataFrame: DataFrame containing plausibility-filtered isoform-unique GRN edges.
    '''
    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_grn = op.join(results_dir, 'grn')

    as_aware_prefix = f"{config['tissue']}_as-aware.network_"
    
    filtered_iso_unique_path = op.join(results_dir_grn, as_aware_prefix + 'plausibility_filtered_iso_unique.tsv')
    grn_filt_iso_unique = pd.read_csv(filtered_iso_unique_path, sep='\t')

    return grn_filt_iso_unique


def do_gsea(config, targets):
    """
    Performs Gene Set Enrichment Analysis (GSEA) on plausibility-filtered isoform-unique GRN targets.

    This function:
    - Loads the isoform-unique GRN edges that passed plausibility filtering.
    - Extracts the top 3000 unique target gene names for enrichment.
    - Filters out NaNs and constructs a background gene set from the expression data.
    - Runs GSEA using gseapy and stores the results.

    Parameters:
        config (dict): Configuration dictionary with paths, including:
            - 'results_dir': Base directory for results.
            - 'tissue': Tissue name used to locate plausibility-filtered GRN data.
        targets (pd.DataFrame): DataFrame with expression data for target genes.
            Must include a column 'gene_id' with Ensembl gene IDs.

    Returns:
        gseapy.enrichr.Enrichr: Enrichment result object containing pathway analysis results.
    """
    #load fully processed grn that underwent plausibility filtering
    iso_unique_f = load_plausibility_filtered(config) 


    # compile target genes for gsea
    target_genes = list(iso_unique_f['target_gene_name'][:7000].unique())
    target_genes = [x for x in target_genes if pd.notnull(x)]
    background_genes = list(targets['gene_id'].drop_duplicates())

    results_dir = op.join(config['results_dir'], config['tissue'])
    results_dir_gsea = op.join(results_dir, 'gsea')

    # do gsea with gseapy
    enr = gp.enrichr(
    gene_list=target_genes,
    gene_sets=[
    'GO_Biological_Process_2021',
    'KEGG_2021_Human',
    'Reactome_2022'
    ],
    background=background_genes,
    organism='human',
    outdir=results_dir_gsea,
    cutoff=0.05
    )

    return enr



def main():
    parser = argparse.ArgumentParser(description="Process a file from the command line.")
    
    # Add the file argument
    parser.add_argument('-f', type=str, help='Config File')
  
    # load config file
    args = parser.parse_args()
    with open(args.f, 'r') as f:
        config = yaml.safe_load(f)

    biomart = pd.read_csv(config['biomart'], sep='\t')
    tf_list = read_tf_list(config['tf_list'], biomart)

    print('Load data')
    transcript_tfs, gene_tfs, targets = magnet.load_magnet_data(config, tf_list)
    #transcript_tfs, gene_tfs, targets = load_data(config, biomart, tf_list)

    do_gsea(config, targets)
    
if __name__ == "__main__":
    main()
    