# Cell 2
import pandas as pd
import numpy as np
import os
import os.path as op
import yaml
# ---

# Cell 3
CONDITION = "DCM"  # Options: 'DCM', 'HCM', 'NF'

# Paths
PART2_RESULTS = "./results_magnet_edge_cate/"
PART3_RESULTS = "./results_magnet_plau_filter/"
BIOMART_PATH = "biomart.txt"

# MAGNet Expression data (for dominance calculation)
MAGNET_TPM_PATH = f"{CONDITION}_magnet_prefiltered_tpm.tsv"

# Annotation resources
APPRIS_PATH = "appris_data.appris.txt"  # APPRIS principal isoform annotations
DIGGER_PATH = "digger_data.csv"  # DIGGER exon/domain annotations

# Filtering Thresholds
DOM_REG = 0.9        # Regulator dominance threshold (one isoform has >90% expression)
DOM_TGT = 0.9        # Target dominance threshold
FC_EQ_MAX = 2.0      # Max allowed fold-change discrepancy for "equivalent" cases
MIN_TX_PER_GENE = 2  # Minimum isoforms for splicing claims
EPSILON = 1e-6

os.makedirs(PART3_RESULTS, exist_ok=True)

# ---

# Cell 5
# Load filtered networks from Part 2
net1 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net1_filtered.tsv"), sep='\t')
net2 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net2_filtered.tsv"), sep='\t')
net3 = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_net3_filtered.tsv"), sep='\t')

# ---

# Cell 6
# Load categorization tables from Part 2 v2
set_a = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_a_source_resolution.tsv"), sep='\t')
set_b = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_b_target_resolution.tsv"), sep='\t')
set_c = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_c_sf_splicing.tsv"), sep='\t')
set_d = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_d_tfsf_joint_ambiguous.tsv"), sep='\t')

# Load unpacked transcript-level tables
set_b_unpacked = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_b_unpacked_txlevel.tsv"), sep='\t')
set_c_unpacked = pd.read_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_c_unpacked_txlevel.tsv"), sep='\t')

# Display Part 2 v2 structure

# ---

# Cell 7
# Load BioMart mappings
biomart = pd.read_csv(BIOMART_PATH, sep='\t')
tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()
gene2symbol = dict(zip(biomart['Gene stable ID'], biomart['Gene name']))

# ---

# Cell 9
# Load MAGNet expression data for dominance calculation
transcript_tpm = pd.read_csv(MAGNET_TPM_PATH, sep='\t')
sample_cols = [c for c in transcript_tpm.columns if c not in ['transcript_id', 'gene_id']]
# ---

# Cell 10
def compute_dominance_metrics(transcript_tpm, sample_cols, min_tpm=0.1):
    """
    Compute dominance metrics for each gene.
    
    Returns:
    - gene_dominance: dict of gene_id -> max isoform expression share
    - gene_n_isoforms: dict of gene_id -> number of expressed isoforms
    - tx_expression_share: dict of tx_id -> expression share within gene
    """
    # Compute mean TPM per transcript
    tx_mean_tpm = transcript_tpm[sample_cols].mean(axis=1)
    transcript_tpm = transcript_tpm.copy()
    transcript_tpm['mean_tpm'] = tx_mean_tpm
    
    # Filter to expressed transcripts
    expressed = transcript_tpm[transcript_tpm['mean_tpm'] >= min_tpm].copy()
    
    # Compute gene totals
    gene_totals = expressed.groupby('gene_id')['mean_tpm'].sum()
    
    # Compute expression share per transcript
    expressed['gene_total'] = expressed['gene_id'].map(gene_totals)
    expressed['expr_share'] = expressed['mean_tpm'] / (expressed['gene_total'] + EPSILON)
    
    # Dominance = max expression share per gene
    gene_dominance = expressed.groupby('gene_id')['expr_share'].max().to_dict()
    
    # Number of expressed isoforms per gene
    gene_n_isoforms = expressed.groupby('gene_id').size().to_dict()
    
    # Transcript expression share
    tx_expression_share = dict(zip(expressed['transcript_id'], expressed['expr_share']))
    
    return gene_dominance, gene_n_isoforms, tx_expression_share
# ---

# Cell 11
gene_dominance, gene_n_isoforms, tx_expr_share = compute_dominance_metrics(
    transcript_tpm, sample_cols
)

# ---

# Cell 14
def filter_set_a(set_a, gene_dominance, gene_n_isoforms,
                  dom_threshold=DOM_REG, fc_eq_max=FC_EQ_MAX, eps=EPSILON):
    """
    Apply plausibility filters to Set A (source resolution).
    Updated for Part 2 v2 column names.
    """
    df = set_a.copy()
    
    # Add dominance and isoform count
    df['reg_dominance'] = df['regulator_gene'].map(gene_dominance).fillna(0.5)
    df['reg_n_isoforms'] = df['regulator_gene'].map(gene_n_isoforms).fillna(1)
    
    # Initialize filter columns
    df['is_plausible'] = True
    df['filter_reasons'] = ''
    
    # Identify importance columns (v2 uses S1_mean/S2_mean for ratios)
    s1_col = 'S1_mean'
    s2_col = 'S2_mean'
    
    # Compute ratio for reference (using mean importance for ratio - AlterNet 1.0 style)
    df['ratio_S2_S1'] = df[s2_col] / (df[s1_col] + eps)
    
    # Filter 1: Regulator equivalence (single isoform)
    mask_single_iso = df['reg_n_isoforms'] == 1
    mask_iso_specific = df['source_category'].isin(['source_isoform_specific', 'source_gene_specific'])
    
    implausible_single = mask_single_iso & mask_iso_specific
    df.loc[implausible_single, 'is_plausible'] = False
    df.loc[implausible_single, 'filter_reasons'] += 'single_isoform_reg;'
    
    # Filter 2: Regulator dominance (one isoform dominates)
    mask_dominant = df['reg_dominance'] >= dom_threshold
    mask_isoform_specific = df['source_category'] == 'source_isoform_specific'
    
    # For isoform-specific with dominant regulator, require stronger evidence
    # (S2 >> S1, ratio should be very high)
    weak_evidence = mask_dominant & mask_isoform_specific & (df['ratio_S2_S1'] < 2.0)
    df.loc[weak_evidence, 'is_plausible'] = False
    df.loc[weak_evidence, 'filter_reasons'] += 'dominant_reg_weak_evidence;'
    
    # Filter 3: Equivalent consistency
    mask_equivalent = df['source_category'] == 'source_equivalent'
    
    # Require ratio within bounds for equivalent
    ratio_ok = (df['ratio_S2_S1'] >= 1/fc_eq_max) & (df['ratio_S2_S1'] <= fc_eq_max)
    ratio_ok = ratio_ok | (df[s1_col] == 0) | (df[s2_col] == 0)  # Allow edge cases
    
    inconsistent_equiv = mask_equivalent & ~ratio_ok
    df.loc[inconsistent_equiv, 'source_category'] = 'source_ambiguous'
    df.loc[inconsistent_equiv, 'filter_reasons'] += 'equiv_ratio_out_of_bounds;'
    
    return df
# ---

# Cell 15
set_a_filtered = filter_set_a(set_a, gene_dominance, gene_n_isoforms)

for cat in ['source_isoform_specific', 'source_gene_specific', 'source_equivalent', 'source_ambiguous']:
    total = (set_a_filtered['source_category'] == cat).sum()
    plausible = ((set_a_filtered['source_category'] == cat) & set_a_filtered['is_plausible']).sum()
# ---

# Cell 17
def filter_set_b(set_b, gene_dominance, gene_n_isoforms,
                  dom_threshold=DOM_TGT, fc_eq_max=FC_EQ_MAX, eps=EPSILON):
    """
    Apply plausibility filters to Set B (target resolution).
    Updated for Part 2 v2 structure with reg_type column.
    """
    df = set_b.copy()
    
    # Add dominance and isoform count for TARGET gene
    df['tgt_dominance'] = df['target_gene'].map(gene_dominance).fillna(0.5)
    df['tgt_n_isoforms'] = df['target_gene'].map(gene_n_isoforms).fillna(1)
    
    # Initialize filter columns
    df['is_plausible'] = True
    df['filter_reasons'] = ''
    
    # Identify importance columns (v2 uses S2_mean/S3_mean_sum for ratios)
    s2_col = 'S2_mean'
    s3_col = 'S3_mean_sum'
    
    # Compute ratio for reference (using mean importance for ratio)
    df['ratio_S3_S2'] = df[s3_col] / (df[s2_col] + eps)
    
    # Filter 1: Target equivalence (single isoform)
    mask_single_iso = df['tgt_n_isoforms'] == 1
    mask_isoform_specific = df['target_category'] == 'target_isoform_specific'
    
    implausible_single = mask_single_iso & mask_isoform_specific
    df.loc[implausible_single, 'is_plausible'] = False
    df.loc[implausible_single, 'filter_reasons'] += 'single_isoform_tgt;'
    
    # Filter 2: Target dominance (one isoform dominates)
    mask_dominant = df['tgt_dominance'] >= dom_threshold
    
    # For isoform-specific with dominant target, the claim is weaker
    # but we allow it if the edge genuinely exists in Net3 but not Net2
    # (E3=True, E2=False is already checked in categorization)
    # Here we flag for review but don't automatically reject
    dominant_isoform_specific = mask_dominant & mask_isoform_specific
    df.loc[dominant_isoform_specific, 'filter_reasons'] += 'dominant_tgt_review;'
    
    # Filter 3: Equivalent consistency
    mask_equivalent = df['target_category'] == 'target_equivalent'
    
    # Require ratio within bounds for equivalent
    # Handle cases where S2 or S3 might be 0
    ratio_ok = (df['ratio_S3_S2'] >= 1/fc_eq_max) & (df['ratio_S3_S2'] <= fc_eq_max)
    ratio_ok = ratio_ok | (df[s2_col] == 0) | (df[s3_col] == 0)  # Allow edge cases
    
    inconsistent_equiv = mask_equivalent & ~ratio_ok
    df.loc[inconsistent_equiv, 'target_category'] = 'target_ambiguous'
    df.loc[inconsistent_equiv, 'filter_reasons'] += 'equiv_ratio_out_of_bounds;'
    
    # Additional: Flag TF_SF (tf_like) edges for annotation
    tfsf_edges = df['reg_type'] == 'TF_SF'
    df.loc[tfsf_edges, 'filter_reasons'] += 'tfsf_tf_like;'
    
    return df
# ---

# Cell 18
set_b_filtered = filter_set_b(set_b, gene_dominance, gene_n_isoforms)

for cat in ['target_isoform_specific', 'target_gene_specific', 'target_equivalent', 'target_ambiguous']:
    total = (set_b_filtered['target_category'] == cat).sum()
    plausible = ((set_b_filtered['target_category'] == cat) & set_b_filtered['is_plausible']).sum()

# Show by reg_type if available
# ---

# Cell 20
def filter_set_c(set_c, gene_n_isoforms, min_tx_per_gene=MIN_TX_PER_GENE):
    """
    Apply plausibility filters to Set C (SF splicing evidence).
    Updated for Part 2 v2 structure with reg_type column.
    """
    df = set_c.copy()
    
    # Add isoform count for target gene
    df['tgt_n_isoforms'] = df['target_gene'].map(gene_n_isoforms).fillna(1)
    
    # Initialize filter columns
    df['is_plausible'] = True
    df['filter_reasons'] = ''
    
    # Filter 1: Multi-isoform requirement for splicing claims
    mask_splicing = df['sf_category'].str.contains('splicing_supported', na=False)
    mask_single_iso = df['tgt_n_isoforms'] < min_tx_per_gene
    
    invalid_splicing = mask_splicing & mask_single_iso
    
    # Reclassify to expression_associated
    df.loc[invalid_splicing, 'sf_category'] = 'sf_expression_associated'
    df.loc[invalid_splicing, 'filter_reasons'] += 'insufficient_isoforms_for_splicing;'
    
    # Filter 2: Verify specific criteria
    # Already enforced in Part 2 categorization, but double-check
    # Specific should have push_pull OR n_sig >= 2
    mask_specific = df['sf_category'] == 'sf_splicing_supported_specific'
        
    has_push_pull = df['push_pull'] == True
        
    has_multi_sig = df['n_sig'] >= 2
        
    # If specific but doesn't meet criteria, downgrade
    invalid_specific = mask_specific & ~has_push_pull & ~has_multi_sig
    df.loc[invalid_specific, 'sf_category'] = 'sf_expression_associated'
    df.loc[invalid_specific, 'filter_reasons'] += 'specific_criteria_not_met;'
    
    # Filter 3: QC enforcement
    unreliable = df['usage_reliable'] == False
    splicing_unreliable = mask_splicing & unreliable
        
    df.loc[splicing_unreliable, 'sf_category'] = 'sf_ambiguous'
    df.loc[splicing_unreliable, 'filter_reasons'] += 'usage_unreliable;'
    
    # Additional: Flag TF_SF (sf_like) edges for annotation
    tfsf_edges = df['reg_type'] == 'TF_SF'
    df.loc[tfsf_edges, 'filter_reasons'] += 'tfsf_sf_like;'
    
    return df
# ---

# Cell 21
set_c_filtered = filter_set_c(set_c, gene_n_isoforms)

for cat in ['sf_splicing_supported_specific', 'sf_expression_associated', 'sf_ambiguous']:
    total = (set_c_filtered['sf_category'] == cat).sum()

# Show by reg_type if available
# ---

# Cell 23
def filter_set_d(set_d, gene_dominance, gene_n_isoforms, 
                  dom_threshold=DOM_TGT, min_tx=MIN_TX_PER_GENE):
    """
    Apply plausibility filters to Set D (TF+SF joint/ambiguous only in v2).
    """
    df = set_d.copy()
    
    # Add target gene info
    df['tgt_dominance'] = df['target_gene'].map(gene_dominance).fillna(0.5)
    df['tgt_n_isoforms'] = df['target_gene'].map(gene_n_isoforms).fillna(1)
    
    # Initialize filter columns
    df['is_plausible'] = True
    df['filter_reasons'] = ''
    
    # Filter 1: QC flag check
    qc_failed = df['qc_ok'] == False
    df.loc[qc_failed, 'tfsf_category'] = 'tfsf_ambiguous'
    df.loc[qc_failed, 'filter_reasons'] += 'qc_failed;'
    
    # Filter 2: Multi-isoform requirement for joint claims
    mask_joint = df['tfsf_category'] == 'tfsf_joint'
    mask_single_iso = df['tgt_n_isoforms'] < min_tx
    
    # Joint claims need multi-isoform targets (SF evidence needs splicing)
    invalid_joint = mask_joint & mask_single_iso
    df.loc[invalid_joint, 'tfsf_category'] = 'tfsf_ambiguous'
    df.loc[invalid_joint, 'filter_reasons'] += 'joint_single_isoform;'
    
    # Filter 3: Target dominance check for joint claims
    mask_dominant_tgt = df['tgt_dominance'] >= dom_threshold
    
    # Flag dominant target joint claims for review
    dominant_joint = mask_dominant_tgt & mask_joint
    df.loc[dominant_joint, 'filter_reasons'] += 'dominant_tgt_joint_claim;'
    
    return df
# ---

# Cell 24
set_d_filtered = filter_set_d(set_d, gene_dominance, gene_n_isoforms)

for cat in ['tfsf_joint', 'tfsf_ambiguous']:
    total = (set_d_filtered['tfsf_category'] == cat).sum()
# ---

# Cell 26
# Propagate plausibility to unpacked tables

# Set B Unpacked

# Build plausibility lookup from Set B: (regulator_tx, target_gene) -> is_plausible
t2_plaus_key = set_b_filtered[['regulator_tx', 'target_gene', 'is_plausible', 'target_category']].copy()
t2_plaus_key['edge_tg'] = t2_plaus_key['regulator_tx'] + '|' + t2_plaus_key['target_gene']
t2_plaus_dict = dict(zip(t2_plaus_key['edge_tg'], t2_plaus_key['is_plausible']))
t2_cat_dict = dict(zip(t2_plaus_key['edge_tg'], t2_plaus_key['target_category']))

set_b_unpacked['edge_tg'] = set_b_unpacked['regulator_tx'] + '|' + set_b_unpacked['target_gene']
set_b_unpacked['is_plausible'] = set_b_unpacked['edge_tg'].map(t2_plaus_dict).fillna(False)
# Update target_category in case Part 3 reclassified some (e.g. equiv -> ambiguous)
set_b_unpacked['target_category'] = set_b_unpacked['edge_tg'].map(t2_cat_dict).fillna(
    set_b_unpacked['target_category'])
set_b_unpacked = set_b_unpacked.drop(columns=['edge_tg'])

# Set C Unpacked

t3_plaus_key = set_c_filtered[['sf_tx', 'target_gene', 'is_plausible', 'sf_category']].copy()
t3_plaus_key['edge_sg'] = t3_plaus_key['sf_tx'] + '|' + t3_plaus_key['target_gene']
t3_plaus_dict = dict(zip(t3_plaus_key['edge_sg'], t3_plaus_key['is_plausible']))
t3_cat_dict = dict(zip(t3_plaus_key['edge_sg'], t3_plaus_key['sf_category']))

set_c_unpacked['edge_sg'] = set_c_unpacked['sf_tx'] + '|' + set_c_unpacked['target_gene']
set_c_unpacked['is_plausible'] = set_c_unpacked['edge_sg'].map(t3_plaus_dict).fillna(False)
set_c_unpacked['sf_category'] = set_c_unpacked['edge_sg'].map(t3_cat_dict).fillna(
    set_c_unpacked['sf_category'])
set_c_unpacked = set_c_unpacked.drop(columns=['edge_sg'])

# ---

# Cell 29
# Load APPRIS annotations
appris = pd.read_csv(APPRIS_PATH, sep='\t')
appris['transcript_id'] = appris['Transcript ID'].str.split('.').str[0]
tx_appris_label = dict(zip(appris['transcript_id'], appris['APPRIS Annotation']))
# ---

# Cell 30
# Load DIGGER annotations
digger = pd.read_csv(DIGGER_PATH, sep=',')
digger['transcript_id'] = digger['Transcript stable ID'].str.split('.').str[0]

# Raw per-transcript lookups
tx_exon_count = digger.groupby('transcript_id')['Exon rank in transcript'].nunique().to_dict()

domain_counts = digger[digger['Pfam ID'].notna()].groupby('transcript_id')['Pfam ID'].nunique()
tx_domain_count = domain_counts.to_dict()

def get_domains(group):
    domains = group['Pfam ID'].dropna().unique()
    return ','.join(sorted(domains)) if len(domains) > 0 else ''
tx_domains = digger.groupby('transcript_id').apply(get_domains).to_dict()

# Set-based lookups (for unique/shared comparison)
tx_exon_set = digger.groupby('transcript_id')['Exon rank in transcript'].apply(
    lambda x: set(x.dropna().unique())
).to_dict()

tx_domain_set = digger[digger['Pfam ID'].notna()].groupby('transcript_id')['Pfam ID'].apply(
    lambda x: set(x.unique())
).to_dict()

# Gene-to-transcript mapping within DIGGER
digger_gene2tx = digger.groupby(
    digger['Gene stable ID'].str.split('.').str[0]
)['transcript_id'].apply(set).to_dict()
# ---

# Cell 31
# Use BioMart for basic transcript annotations

# Gene type (BioMart has 'Gene type', not transcript-level biotype)
gene_biotype = dict(zip(biomart['Gene stable ID'], biomart['Gene type']))
tx_biotype = {tx: gene_biotype.get(gene, '') for tx, gene in tx2gene.items()}

# Gene symbol
tx_gene_symbol = {tx: gene2symbol.get(gene, '')
                  for tx, gene in tx2gene.items()}
# ---

# Cell 32
# DIGGER comparison functions (unique / shared features)

def get_principal_tx(gene_id):
    """Return the principal isoform for a gene using APPRIS labels."""
    candidates = digger_gene2tx.get(gene_id, set())
    # Among transcripts of this gene that have APPRIS labels,
    # find the one labeled 'PRINCIPAL:1' (or the highest-ranked principal)
    best_tx, best_rank = None, 999
    for tx in candidates:
        label = tx_appris_label.get(tx, '')
        if 'PRINCIPAL' in str(label).upper():
            # APPRIS ranks: PRINCIPAL:1 > PRINCIPAL:2 > ... > ALTERNATIVE:1 > ...
            parts = str(label).split(':')
            rank = int(parts[-1]) if len(parts) > 1 and parts[-1].isdigit() else 99
            if rank < best_rank:
                best_rank = rank
                best_tx = tx
    return best_tx

def get_unique_features(tx_id, gene_id):
    """
    Compare tx_id against a baseline and return unique exons/domains.
    Baseline = principal isoform (if tx_id is not principal),
               or all other isoforms (if tx_id is principal or no principal exists).
    Returns: (unique_exons_str, unique_domains_str)
    """
    if pd.isna(tx_id) or tx_id == '':
        return '', ''

    tx_exons = tx_exon_set.get(tx_id, set())
    tx_doms  = tx_domain_set.get(tx_id, set())
    if not tx_exons and not tx_doms:
        return '', ''

    gene_txs = digger_gene2tx.get(gene_id, set()) if gene_id else set()
    principal = get_principal_tx(gene_id)

    if principal and principal != tx_id:
        # Compare against principal
        baseline_exons = tx_exon_set.get(principal, set())
        baseline_doms  = tx_domain_set.get(principal, set())
    else:
        # tx_id IS principal (or no principal) -> compare against union of all others
        others = gene_txs - {tx_id}
        baseline_exons = set().union(*(tx_exon_set.get(t, set()) for t in others)) if others else set()
        baseline_doms  = set().union(*(tx_domain_set.get(t, set()) for t in others)) if others else set()

    unique_exons  = tx_exons - baseline_exons
    unique_doms   = tx_doms  - baseline_doms

    ue_str = ','.join(str(e) for e in sorted(unique_exons)) if unique_exons else ''
    ud_str = ','.join(sorted(unique_doms)) if unique_doms else ''
    return ue_str, ud_str

def get_shared_features(gene_id):
    """
    Return exons and domains shared by ALL isoforms of a gene (intersection).
    Returns: (shared_exons_str, shared_domains_str)
    """
    if pd.isna(gene_id) or gene_id == '':
        return '', ''

    gene_txs = digger_gene2tx.get(gene_id, set())
    if not gene_txs:
        return '', ''

    # Only consider transcripts that appear in DIGGER
    txs_with_data = [t for t in gene_txs if t in tx_exon_set]
    if not txs_with_data:
        return '', ''

    shared_exons = tx_exon_set.get(txs_with_data[0], set()).copy()
    shared_doms  = tx_domain_set.get(txs_with_data[0], set()).copy()
    for t in txs_with_data[1:]:
        shared_exons &= tx_exon_set.get(t, set())
        shared_doms  &= tx_domain_set.get(t, set())

    se_str = ','.join(str(e) for e in sorted(shared_exons)) if shared_exons else ''
    sd_str = ','.join(sorted(shared_doms)) if shared_doms else ''
    return se_str, sd_str

def annotate_digger_raw(tx_id):
    """Return basic DIGGER annotation for a single transcript."""
    if pd.isna(tx_id) or tx_id == '':
        return '', '', ''
    return (
        tx_exon_count.get(tx_id, ''),
        tx_domain_count.get(tx_id, ''),
        tx_domains.get(tx_id, '')
    )
# ---

# Cell 34
# Annotate Set A (regulator annotations via best_tx)

set_a_filtered['reg_appris'] = set_a_filtered['best_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_a_filtered['reg_biotype'] = set_a_filtered['best_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)
set_a_filtered['reg_expr_share'] = set_a_filtered['best_tx'].map(
    lambda x: tx_expr_share.get(x, np.nan) if pd.notna(x) else np.nan
)

# Add gene symbols
set_a_filtered['reg_symbol'] = set_a_filtered['regulator_gene'].map(gene2symbol)
set_a_filtered['tgt_symbol'] = set_a_filtered['target_gene'].map(gene2symbol)

# DIGGER annotation for regulator (best_tx)
    raw = set_a_filtered['best_tx'].apply(annotate_digger_raw)
    set_a_filtered['reg_exon_count']  = raw.apply(lambda x: x[0])
    set_a_filtered['reg_domain_count'] = raw.apply(lambda x: x[1])
    set_a_filtered['reg_domains']     = raw.apply(lambda x: x[2])

    # Category-dependent: unique for source_isoform_specific, shared for source_gene_specific
    def digger_set_a(row):
        cat = row.get('source_category', '')
        tx  = row.get('best_tx', '')
        gene = row.get('regulator_gene', '')
        if cat == 'source_isoform_specific':
            ue, ud = get_unique_features(tx, gene)
            return pd.Series({'reg_unique_exons': ue, 'reg_unique_domains': ud,
                              'reg_shared_exons': '', 'reg_shared_domains': ''})
        elif cat == 'source_gene_specific':
            se, sd = get_shared_features(gene)
            return pd.Series({'reg_unique_exons': '', 'reg_unique_domains': '',
                              'reg_shared_exons': se, 'reg_shared_domains': sd})
        return pd.Series({'reg_unique_exons': '', 'reg_unique_domains': '',
                          'reg_shared_exons': '', 'reg_shared_domains': ''})

    digger_cols = set_a_filtered.apply(digger_set_a, axis=1)
    set_a_filtered = pd.concat([set_a_filtered, digger_cols], axis=1)
# ---

# Cell 35
# Annotate Set B (target annotations)

# Regulator transcript annotations
set_b_filtered['reg_appris'] = set_b_filtered['regulator_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_b_filtered['reg_biotype'] = set_b_filtered['regulator_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Gene symbols
set_b_filtered['reg_symbol'] = set_b_filtered['regulator_gene'].map(gene2symbol)
set_b_filtered['tgt_symbol'] = set_b_filtered['target_gene'].map(gene2symbol)

# DIGGER: raw annotation for regulator
raw = set_b_filtered['regulator_tx'].apply(annotate_digger_raw)
set_b_filtered['reg_exon_count']  = raw.apply(lambda x: x[0])
set_b_filtered['reg_domain_count'] = raw.apply(lambda x: x[1])
set_b_filtered['reg_domains']     = raw.apply(lambda x: x[2])

# Category-dependent for target side (using target_tx_resolved)
def digger_set_b(row):
    cat = row.get('target_category', '')
    tx  = row.get('target_tx_resolved', '')
    gene = row.get('target_gene', '')
    if cat == 'target_isoform_specific':
        ue, ud = get_unique_features(tx, gene)
        return pd.Series({'tgt_unique_exons': ue, 'tgt_unique_domains': ud,
                          'tgt_shared_exons': '', 'tgt_shared_domains': ''})
    elif cat == 'target_gene_specific':
        se, sd = get_shared_features(gene)
        return pd.Series({'tgt_unique_exons': '', 'tgt_unique_domains': '',
                          'tgt_shared_exons': se, 'tgt_shared_domains': sd})
    return pd.Series({'tgt_unique_exons': '', 'tgt_unique_domains': '',
                      'tgt_shared_exons': '', 'tgt_shared_domains': ''})

digger_cols = set_b_filtered.apply(digger_set_b, axis=1)
set_b_filtered = pd.concat([set_b_filtered, digger_cols], axis=1)
# ---

# Cell 36
# Annotate Set C (SF and target annotations)

# SF transcript annotations
set_c_filtered['sf_appris'] = set_c_filtered['sf_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_c_filtered['sf_biotype'] = set_c_filtered['sf_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Gene symbols
set_c_filtered['sf_symbol'] = set_c_filtered['sf_gene'].map(gene2symbol)
set_c_filtered['tgt_symbol'] = set_c_filtered['target_gene'].map(gene2symbol)

# DIGGER: raw annotation for SF transcript
raw = set_c_filtered['sf_tx'].apply(annotate_digger_raw)
set_c_filtered['sf_exon_count']  = raw.apply(lambda x: x[0])
set_c_filtered['sf_domain_count'] = raw.apply(lambda x: x[1])
set_c_filtered['sf_domains']     = raw.apply(lambda x: x[2])

# Category-dependent for target side (using best_tx from splicing correlation)
def digger_set_c(row):
    cat = row.get('sf_category', '')
    tx  = row.get('best_tx', '') if pd.notna(row.get('best_tx', '')) else row.get('target_tx_resolved', '')
    gene = row.get('target_gene', '')
    if cat == 'sf_splicing_supported_specific':
        ue, ud = get_unique_features(tx, gene)
        return pd.Series({'tgt_unique_exons': ue, 'tgt_unique_domains': ud,
                          'tgt_shared_exons': '', 'tgt_shared_domains': ''})
    elif cat == 'sf_splicing_supported_broad':
        se, sd = get_shared_features(gene)
        return pd.Series({'tgt_unique_exons': '', 'tgt_unique_domains': '',
                          'tgt_shared_exons': se, 'tgt_shared_domains': sd})
    return pd.Series({'tgt_unique_exons': '', 'tgt_unique_domains': '',
                      'tgt_shared_exons': '', 'tgt_shared_domains': ''})

digger_cols = set_c_filtered.apply(digger_set_c, axis=1)
set_c_filtered = pd.concat([set_c_filtered, digger_cols], axis=1)
# ---

# Cell 37
# Annotate Set D (both regulator and target annotations)

# Regulator annotations
set_d_filtered['reg_appris'] = set_d_filtered['reg_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_d_filtered['reg_biotype'] = set_d_filtered['reg_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Target annotations
set_d_filtered['tgt_appris'] = set_d_filtered['target_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_d_filtered['tgt_biotype'] = set_d_filtered['target_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Gene symbols
set_d_filtered['reg_symbol'] = set_d_filtered['reg_gene'].map(gene2symbol)
set_d_filtered['tgt_symbol'] = set_d_filtered['target_gene'].map(gene2symbol)

# DIGGER: raw annotation for unpacked tables
# Set B unpacked
raw = set_b_unpacked['regulator_tx'].apply(annotate_digger_raw)
set_b_unpacked['reg_exon_count']  = raw.apply(lambda x: x[0])
set_b_unpacked['reg_domain_count'] = raw.apply(lambda x: x[1])
set_b_unpacked['reg_domains']     = raw.apply(lambda x: x[2])

raw = set_b_unpacked['target_tx'].apply(annotate_digger_raw)
set_b_unpacked['tgt_exon_count']  = raw.apply(lambda x: x[0])
set_b_unpacked['tgt_domain_count'] = raw.apply(lambda x: x[1])
set_b_unpacked['tgt_domains']     = raw.apply(lambda x: x[2])

# Set C unpacked
raw = set_c_unpacked['sf_tx'].apply(annotate_digger_raw)
set_c_unpacked['sf_exon_count']  = raw.apply(lambda x: x[0])
set_c_unpacked['sf_domain_count'] = raw.apply(lambda x: x[1])
set_c_unpacked['sf_domains']     = raw.apply(lambda x: x[2])

raw = set_c_unpacked['target_tx'].apply(annotate_digger_raw)
set_c_unpacked['tgt_exon_count']  = raw.apply(lambda x: x[0])
set_c_unpacked['tgt_domain_count'] = raw.apply(lambda x: x[1])
set_c_unpacked['tgt_domains']     = raw.apply(lambda x: x[2])
# ---

# Cell 39
# Annotate unpacked tables

# Set B Unpacked

# Regulator annotations
set_b_unpacked['reg_appris'] = set_b_unpacked['regulator_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_b_unpacked['reg_biotype'] = set_b_unpacked['regulator_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Target transcript annotations
set_b_unpacked['tgt_appris'] = set_b_unpacked['target_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_b_unpacked['tgt_biotype'] = set_b_unpacked['target_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Gene symbols
set_b_unpacked['reg_symbol'] = set_b_unpacked['regulator_gene'].map(gene2symbol)
set_b_unpacked['tgt_symbol'] = set_b_unpacked['target_gene'].map(gene2symbol)

# Set C Unpacked

# SF annotations
set_c_unpacked['sf_appris'] = set_c_unpacked['sf_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_c_unpacked['sf_biotype'] = set_c_unpacked['sf_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Target transcript annotations
set_c_unpacked['tgt_appris'] = set_c_unpacked['target_tx'].map(
    lambda x: tx_appris_label.get(x, '') if pd.notna(x) else ''
)
set_c_unpacked['tgt_biotype'] = set_c_unpacked['target_tx'].map(
    lambda x: tx_biotype.get(x, '') if pd.notna(x) else ''
)

# Gene symbols
set_c_unpacked['sf_symbol'] = set_c_unpacked['sf_gene'].map(gene2symbol)
set_c_unpacked['tgt_symbol'] = set_c_unpacked['target_gene'].map(gene2symbol)

# Also annotate presentation columns on Set B and Set C
for col_name in ['target_tx_resolved', 'target_tx_dominant']:
    set_b_filtered[f'{col_name}_appris'] = set_b_filtered[col_name].map(
        lambda x: tx_appris_label.get(x, '') if pd.notna(x) and x != '' else ''
    )
    set_b_filtered[f'{col_name}_biotype'] = set_b_filtered[col_name].map(
        lambda x: tx_biotype.get(x, '') if pd.notna(x) and x != '' else ''
    )

for col_name in ['target_tx_resolved', 'target_tx_dominant']:
    set_c_filtered[f'{col_name}_appris'] = set_c_filtered[col_name].map(
        lambda x: tx_appris_label.get(x, '') if pd.notna(x) and x != '' else ''
    )
    set_c_filtered[f'{col_name}_biotype'] = set_c_filtered[col_name].map(
        lambda x: tx_biotype.get(x, '') if pd.notna(x) and x != '' else ''
    )

# ---

# Cell 41
# Save filtered and annotated tables

# Save all rows (with plausibility flags)
set_a_filtered.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_a_filtered_annotated.tsv"), 
                       sep='\t', index=False)
set_b_filtered.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_b_filtered_annotated.tsv"), 
                       sep='\t', index=False)
set_c_filtered.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_c_filtered_annotated.tsv"), 
                       sep='\t', index=False)
set_d_filtered.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_d_filtered_annotated.tsv"), 
                       sep='\t', index=False)

# Save unpacked tables (all rows with plausibility flags)
set_b_unpacked.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_b_unpacked_filtered_annotated.tsv"),
                       sep='\t', index=False)
set_c_unpacked.to_csv(op.join(PART3_RESULTS, f"{CONDITION}_set_c_unpacked_filtered_annotated.tsv"),
                       sep='\t', index=False)

# Save plausible-only versions
set_a_filtered[set_a_filtered['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_a_plausible.tsv"), sep='\t', index=False)
set_b_filtered[set_b_filtered['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_b_plausible.tsv"), sep='\t', index=False)
set_c_filtered[set_c_filtered['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_c_plausible.tsv"), sep='\t', index=False)
set_d_filtered[set_d_filtered['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_d_plausible.tsv"), sep='\t', index=False)

# Save plausible unpacked tables
set_b_unpacked[set_b_unpacked['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_b_unpacked_plausible.tsv"), sep='\t', index=False)
set_c_unpacked[set_c_unpacked['is_plausible']].to_csv(
    op.join(PART3_RESULTS, f"{CONDITION}_set_c_unpacked_plausible.tsv"), sep='\t', index=False)

# ---

# Cell 42
# Save summary YAML
summary = {
    'condition': CONDITION,
    'version': 'v2',
    'parameters': {
        'DOM_REG': DOM_REG,
        'DOM_TGT': DOM_TGT,
        'FC_EQ_MAX': FC_EQ_MAX,
        'MIN_TX_PER_GENE': MIN_TX_PER_GENE
    },
    'set_a': {
        'total': len(set_a_filtered),
        'plausible': int(set_a_filtered['is_plausible'].sum()),
        'categories': {k: int(v) for k, v in set_a_filtered[set_a_filtered['is_plausible']]['source_category'].value_counts().to_dict().items()}
    },
    'set_b': {
        'total': len(set_b_filtered),
        'plausible': int(set_b_filtered['is_plausible'].sum()),
        'tf_edges': int((set_b_filtered['reg_type'] == 'TF').sum()),
        'tfsf_tf_like_edges': int((set_b_filtered['reg_type'] == 'TF_SF').sum()),
        'categories': {k: int(v) for k, v in set_b_filtered[set_b_filtered['is_plausible']]['target_category'].value_counts().to_dict().items()}
    },
    'set_b_unpacked': {
        'total': int(len(set_b_unpacked)),
        'plausible': int(set_b_unpacked['is_plausible'].sum()),
    },
    'set_c': {
        'total': len(set_c_filtered),
        'plausible': int(set_c_filtered['is_plausible'].sum()),
        'sf_edges': int((set_c_filtered['reg_type'] == 'SF').sum()),
        'tfsf_sf_like_edges': int((set_c_filtered['reg_type'] == 'TF_SF').sum()),
        'categories': {k: int(v) for k, v in set_c_filtered[set_c_filtered['is_plausible']]['sf_category'].value_counts().to_dict().items()}
    },
    'set_c_unpacked': {
        'total': int(len(set_c_unpacked)),
        'plausible': int(set_c_unpacked['is_plausible'].sum()),
    },
    'set_d': {
        'total': len(set_d_filtered),
        'plausible': int(set_d_filtered['is_plausible'].sum()),
        'note': 'Only contains tfsf_joint and tfsf_ambiguous in v2',
        'categories': {k: int(v) for k, v in set_d_filtered[set_d_filtered['is_plausible']]['tfsf_category'].value_counts().to_dict().items()}
    },
    'annotation_sources': {
        'appris': True,
        'digger': True,
        'expression': True
    }
}

with open(op.join(PART3_RESULTS, f"{CONDITION}_part2p5_summary_v2.yaml"), 'w') as f:
    yaml.dump(summary, f, default_flow_style=False)

# ---