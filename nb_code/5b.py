# Cell 2
import pandas as pd
import numpy as np
import os
import os.path as op
import yaml
from collections import defaultdict
# ---

# Cell 3
CONDITION = "NF"  # Options: 'DCM', 'HCM', 'NF'

# Paths
PART3_RESULTS = "./results_magnet_plau_filter/"
PART2_RESULTS = "./results_magnet_edge_cate/"
GPROFILER_RESULTS = "./results_magnet_gprofiler_input/"
BIOMART_PATH = "biomart.txt"

# ORA parameters
TOP_K = 500
PROJECTION_METHOD = 'max'
USE_ORDERED_QUERY = True
ENRICHMENT_PVAL = 0.05

os.makedirs(GPROFILER_RESULTS, exist_ok=True)
# ---

# Cell 5
# Load curated tables from Part 3

set_a = pd.read_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_a_plausible.tsv"), sep='\t')
set_b = pd.read_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_b_plausible.tsv"), sep='\t')

# ---

# Cell 6
# Load filtered networks from Part 2

net1 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net1_filtered.tsv"), sep='\t')
net2 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net2_filtered.tsv"), sep='\t')
net3 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net3_filtered.tsv"), sep='\t')

# ---

# Cell 7
# Load BioMart mappings

biomart = pd.read_csv(BIOMART_PATH, sep='\t')
tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
gene2symbol = dict(zip(biomart['Gene stable ID'], biomart['Gene name']))
symbol2gene = dict(zip(biomart['Gene name'], biomart['Gene stable ID']))

# ---

# Cell 9
def score_targets(edges, target_col, importance_col='median_importance'):
    """
    Compute target scores as weighted in-degree.
    score(target) = sum of importance for all edges targeting it
    """
    scores = edges.groupby(target_col)[importance_col].sum()
    return scores.sort_values(ascending=False)
# ---

# Cell 10
def project_tx_to_gene(tx_scores, tx2gene, method='max'):
    """
    Project transcript-level scores to gene-level.
    
    Parameters:
    - tx_scores: pd.Series mapping transcript_id -> score
    - tx2gene: dict mapping transcript_id -> gene_id
    - method: 'max' (default) or 'sum'
    
    Returns:
    - gene_scores: pd.Series mapping gene_id -> score
    - rep_tx: dict mapping gene_id -> representative transcript
    """
    df = pd.DataFrame({
        'transcript_id': tx_scores.index,
        'score': tx_scores.values
    })
    df['gene_id'] = df['transcript_id'].map(tx2gene)
    df = df.dropna(subset=['gene_id'])
    
    if method == 'max':
        idx = df.groupby('gene_id')['score'].idxmax()
        result = df.loc[idx].set_index('gene_id')
        gene_scores = result['score'].sort_values(ascending=False)
        rep_tx = result['transcript_id'].to_dict()
    elif method == 'sum':
        gene_scores = df.groupby('gene_id')['score'].sum().sort_values(ascending=False)
        rep_tx = {}
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return gene_scores, rep_tx
# ---

# Cell 11
def build_target_list(edges, target_col, importance_col='median_importance',
                      target_type='gene', tx2gene=None, gene2symbol=None,
                      K=500, projection_method='max'):
    """
    Build a ranked target list from an edge set.
    
    Parameters:
    - edges: DataFrame with edges
    - target_col: column containing target IDs
    - importance_col: column containing importance scores
    - target_type: 'gene' or 'transcript'
    - tx2gene: transcript to gene mapping (required if target_type='transcript')
    - gene2symbol: gene ID to symbol mapping
    - K: number of top genes to return
    - projection_method: 'max' or 'sum' (for transcript targets)
    
    Returns:
    - top_df: DataFrame with ranked targets
    - gene_symbols: list of gene symbols for g:Profiler
    - metadata: dict with scoring info
    """
    if len(edges) == 0:
        return pd.DataFrame(), [], {'n_edges': 0}
    
    # Score targets
    target_scores = score_targets(edges, target_col, importance_col)
    
    # Project if transcript targets
    if target_type == 'transcript':
        if tx2gene is None:
            raise ValueError("tx2gene mapping required for transcript targets")
        gene_scores, rep_tx = project_tx_to_gene(target_scores, tx2gene, method=projection_method)
    else:
        gene_scores = target_scores
        rep_tx = {}
    
    # Get top K
    top_genes = gene_scores.head(K)
    
    # Build result DataFrame
    top_df = pd.DataFrame({
        'rank': range(1, len(top_genes) + 1),
        'gene_id': top_genes.index,
        'score': top_genes.values
    })
    
    if gene2symbol is not None:
        top_df['gene_symbol'] = top_df['gene_id'].map(gene2symbol)
        top_df = top_df.dropna(subset=['gene_symbol'])
        gene_symbols = top_df['gene_symbol'].tolist()
    else:
        gene_symbols = top_df['gene_id'].tolist()
    
    if len(rep_tx) > 0:
        top_df['rep_transcript'] = top_df['gene_id'].map(rep_tx)
    
    metadata = {
        'n_edges': len(edges),
        'n_targets': len(target_scores),
        'n_genes': len(gene_scores),
        'n_top_symbols': len(gene_symbols),
        'projection_method': projection_method if target_type == 'transcript' else 'none'
    }
    
    return top_df, gene_symbols, metadata
# ---

# Cell 13
l1_top, l1_symbols, l1_meta = build_target_list(
    edges=net1,
    target_col='target_gene',
    importance_col='median_importance',
    target_type='gene',
    gene2symbol=gene2symbol,
    K=TOP_K
)

# ---

# Cell 15
# Filter Set A for source_isoform_specific edges
l2_edges = set_a[set_a['source_category'] == 'source_isoform_specific'].copy()

importance_col_l2 = 'S2_median'

l2_top, l2_symbols, l2_meta = build_target_list(
    edges=l2_edges,
    target_col='target_gene',
    importance_col=importance_col_l2,
    target_type='gene',
    gene2symbol=gene2symbol,
    K=TOP_K
)
# ---

# Cell 17
# Filter Set B for target_isoform_unique edges
l3_table2 = set_b[set_b['target_category'] == 'target_isoform_specific'].copy()

# Get corresponding Net3 edges for transcript-level scoring
# Set B is at (regulator_tx, target_gene) level
# We need to link to Net3 which has (source_transcript, target_transcript)

# Strategy: Filter Net3 to TF-like edges and get transcript targets
net3_tf_like = net3[net3['reg_type'].isin(['TF', 'TF_SF'])].copy()

# Create lookup from Set B
l3_pairs = set(zip(l3_table2['regulator_tx'], l3_table2['target_gene']))

# Add target_gene to Net3 for matching
net3_tf_like['target_gene'] = net3_tf_like['target_transcript'].map(tx2gene)
net3_tf_like['pair'] = list(zip(net3_tf_like['source_transcript'], net3_tf_like['target_gene']))

# Filter to edges matching Set B target_isoform_unique
l3_edges = net3_tf_like[net3_tf_like['pair'].isin(l3_pairs)].copy()

# Build target list - targets are transcripts, project to genes
l3_top, l3_symbols, l3_meta = build_target_list(
    edges=l3_edges,
    target_col='target_transcript',
    importance_col='median_importance',
    target_type='transcript',
    tx2gene=tx2gene,
    gene2symbol=gene2symbol,
    K=TOP_K,
    projection_method=PROJECTION_METHOD
)

# ---

# Cell 19
# Compare overlap between lists
set_l1 = set(l1_symbols)
set_l2 = set(l2_symbols)
set_l3 = set(l3_symbols)

all_three = set_l1 & set_l2 & set_l3

# ---

# Cell 21
# Save target lists as TSV (with scores)
l1_top.to_csv(op.join(GPROFILER_RESULTS, f"{CONDITION}_L1_canonical.tsv"), sep='\t', index=False)
l2_top.to_csv(op.join(GPROFILER_RESULTS, f"{CONDITION}_L2_source_isoform.tsv"), sep='\t', index=False)
l3_top.to_csv(op.join(GPROFILER_RESULTS, f"{CONDITION}_L3_target_isoform.tsv"), sep='\t', index=False)
# ---

# Cell 22
# Save simple gene lists for web paste (newline-separated)
def save_gene_list(symbols, filepath):
    """Save gene list as newline-separated text file."""
    with open(filepath, 'w') as f:
        f.write('\n'.join(symbols))
    return len(symbols)

n1 = save_gene_list(l1_symbols, op.join(GPROFILER_RESULTS, f"{CONDITION}_L1_genes.txt"))
n2 = save_gene_list(l2_symbols, op.join(GPROFILER_RESULTS, f"{CONDITION}_L2_genes.txt"))
n3 = save_gene_list(l3_symbols, op.join(GPROFILER_RESULTS, f"{CONDITION}_L3_genes.txt"))

# ---

# Cell 23
# Create combined file for multi-query comparison
# Format: Each query on a new line, prefixed with ">"

combined_path = op.join(GPROFILER_RESULTS, f"{CONDITION}_all_lists_combined.txt")

with open(combined_path, 'w') as f:
    f.write(f">L1_Canonical\n")
    f.write(' '.join(l1_symbols) + '\n')
    f.write(f">L2_Source_Isoform\n")
    f.write(' '.join(l2_symbols) + '\n')
    f.write(f">L3_Target_Isoform\n")
    f.write(' '.join(l3_symbols) + '\n')

# ---