# Cell 2
import pandas as pd
import numpy as np
import os
import os.path as op
import yaml
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
# ---

# Cell 3
CONDITION = "DCM"  # Options: 'DCM', 'HCM', 'NF'

# Paths
PART1_RESULTS = "./results_magnet_net_infer/"
PART2_RESULTS = "./results_magnet_edge_cate/"
BIOMART_PATH = "biomart.txt"

# MAGNet Expression data (pre-filtered TPM - same file used in Part 1)
MAGNET_TPM_PATH = f"{CONDITION}_magnet_prefiltered_tpm.tsv"

# Filtering Parameters
MIN_FREQUENCY = 10
IMPORTANCE_PERCENTILE = 0.7  # Keep top 30%

# Categorization Thresholds
R_ISO = 1.5    # Fold change threshold for isoform-specific
R_GENE = 1.5   # Fold change threshold for gene-specific
R_EQ = 1.2     # Equivalence band
EPSILON = 1e-6

# Set C and Set D: PSI/Usage thresholds
GENE_TPM_MIN = 1.0
RHO_MIN = 0.3
Q_MIN = 0.05
DU_MIN = 0.1
DOM_MIN = 0.5
DOM_EQ_MIN = 0.7

# Set C specific
MIN_ISOFORMS_FOR_SPLICING = 2
TOP_M_EXPRESSED = 3

os.makedirs(PART2_RESULTS, exist_ok=True)

# ---

# Cell 5
# Load three networks from Part 1
net1_raw = pd.read_csv(op.join(PART1_RESULTS, f"{CONDITION}_canonical_raw.tsv"), sep='\t')
net2_raw = pd.read_csv(op.join(PART1_RESULTS, f"{CONDITION}_as_aware_source_raw.tsv"), sep='\t')
net3_raw = pd.read_csv(op.join(PART1_RESULTS, f"{CONDITION}_fully_as_aware_raw.tsv"), sep='\t')

print(f"Network 1 (Canonical):      {len(net1_raw):>12,} edges")
print(f"Network 2 (AS-aware src):   {len(net2_raw):>12,} edges")
print(f"Network 3 (Fully AS-aware): {len(net3_raw):>12,} edges")
# ---

# Cell 6
# Load regulator classification
regulator_list = pd.read_csv(op.join(PART1_RESULTS, f"{CONDITION}_regulator_list.csv"))
print(f"\nRegulator classification:")
print(regulator_list['Regulator_type'].value_counts())

tx_to_regtype = dict(zip(regulator_list['Transcript stable ID'], regulator_list['Regulator_type']))
gene_to_regtype = regulator_list.groupby('Gene stable ID')['Regulator_type'].first().to_dict()
# ---

# Cell 7
# Load BioMart
biomart = pd.read_csv(BIOMART_PATH, sep='\t')
tx2gene = dict(zip(biomart['Transcript stable ID'], biomart['Gene stable ID']))
gene2tx = biomart.groupby('Gene stable ID')['Transcript stable ID'].apply(set).to_dict()
print(f"BioMart: {len(tx2gene):,} transcripts, {len(gene2tx):,} genes")
# ---

# Cell 9
net1_raw = net1_raw.rename(columns={'source': 'source_gene', 'target': 'target_gene'})

net2_raw = net2_raw.rename(columns={'source': 'source_transcript', 'target': 'target_gene'})
net2_raw['source_gene'] = net2_raw['source_transcript'].map(tx2gene)

net3_raw = net3_raw.rename(columns={'source': 'source_transcript', 'target': 'target_transcript'})
net3_raw['source_gene'] = net3_raw['source_transcript'].map(tx2gene)
net3_raw['target_gene'] = net3_raw['target_transcript'].map(tx2gene)
# ---

# Cell 11
def filter_edges(df, min_frequency=10, importance_percentile=0.7,
                 freq_col='frequency', imp_col='median_importance'):
    """
    Filter edges by frequency and importance (AlterNet 1.0 style).
    Uses MEDIAN importance for filtering (more robust to outliers).
    """
    n_before = len(df)
    
    # Frequency filter first
    freq_mask = df[freq_col] >= min_frequency
    df_filtered = df[freq_mask].copy()
    n_after_freq = len(df_filtered)
    
    # Importance filter (top percentile) - using median_importance
    importance_threshold = df_filtered[imp_col].quantile(importance_percentile)
    df_filtered = df_filtered[df_filtered[imp_col] >= importance_threshold].copy()
    n_after_imp = len(df_filtered)
    
    print(f"  Before: {n_before:,}")
    print(f"  After frequency filter (>={min_frequency}): {n_after_freq:,}")
    print(f"  After importance filter (top {(1-importance_percentile)*100:.0f}%, threshold={importance_threshold:.4f}): {n_after_imp:,}")
    
    return df_filtered
# ---

# Cell 12
print("Network 1 (Canonical)")
net1_filtered = filter_edges(net1_raw, MIN_FREQUENCY, IMPORTANCE_PERCENTILE)

print("\nNetwork 2 (AS-Aware Source)")
net2_filtered = filter_edges(net2_raw, MIN_FREQUENCY, IMPORTANCE_PERCENTILE)

print("\nNetwork 3 (Fully AS-Aware)")
net3_filtered = filter_edges(net3_raw, MIN_FREQUENCY, IMPORTANCE_PERCENTILE)
# ---

# Cell 14
net1_filtered['reg_type'] = net1_filtered['source_gene'].map(gene_to_regtype)
net2_filtered['reg_type'] = net2_filtered['source_transcript'].map(tx_to_regtype)
net3_filtered['reg_type'] = net3_filtered['source_transcript'].map(tx_to_regtype)

# Split Network 3 by regulator type
net3_tf = net3_filtered[net3_filtered['reg_type'] == 'TF'].copy()
net3_sf = net3_filtered[net3_filtered['reg_type'] == 'SF'].copy()
net3_tfsf = net3_filtered[net3_filtered['reg_type'] == 'TF_SF'].copy()

print(f"Network 3 split by regulator type:")
print(f"  TF edges:    {len(net3_tf):>10,}")
print(f"  SF edges:    {len(net3_sf):>10,}")
print(f"  TF_SF edges: {len(net3_tfsf):>10,}")
# ---

# Cell 16
# Load MAGNet expression data (already condition-specific)
transcript_tpm = pd.read_csv(MAGNET_TPM_PATH, sep='\t')

# Get sample columns (everything except transcript_id and gene_id)
id_cols = {'gene_id', 'transcript_id'}
sample_cols = [c for c in transcript_tpm.columns 
               if c not in id_cols and transcript_tpm[c].dtype in ['float64', 'float32', 'int64', 'int32']]

tpm_numeric = transcript_tpm[sample_cols].astype(float)
tpm_with_gene = tpm_numeric.copy()
tpm_with_gene['_gene_id_'] = transcript_tpm['gene_id'].values

gene_totals = tpm_with_gene.groupby('_gene_id_')[sample_cols].transform('sum')
usage_values = tpm_numeric / (gene_totals + EPSILON)

usage_df = transcript_tpm[['transcript_id', 'gene_id']].copy()
usage_df[sample_cols] = usage_values
reliability_df = (gene_totals >= GENE_TPM_MIN).astype(float)

usage_df_indexed = usage_df.set_index('transcript_id')[sample_cols]
reliability_df_indexed = reliability_df.copy()
reliability_df_indexed.index = usage_df['transcript_id'].values
reliability_df_indexed = reliability_df_indexed[sample_cols]

gene_to_all_transcripts = usage_df.groupby('gene_id')['transcript_id'].apply(list).to_dict()
transcript_mean_expr = transcript_tpm.set_index('transcript_id')[sample_cols].mean(axis=1).to_dict()
expr_df_indexed = transcript_tpm.set_index('transcript_id')[sample_cols]
# ---

# Cell 19
# Filter to TF-like regulators for Set A
net1_tf = net1_filtered[net1_filtered['reg_type'].isin(['TF', 'TF_SF'])].copy()
net2_tf = net2_filtered[net2_filtered['reg_type'].isin(['TF', 'TF_SF'])].copy()

print(f"TF-like edges for Set A: Net1={len(net1_tf):,}, Net2={len(net2_tf):,}")
# ---

# Cell 20
net1_tf['edge_gg'] = net1_tf['source_gene'] + '|' + net1_tf['target_gene']
net2_tf['edge_gg'] = net2_tf['source_gene'] + '|' + net2_tf['target_gene']

# Use MEAN importance for ratio calculation, MEDIAN for ranking
# Net1: aggregate to gene-gene level
net1_agg = net1_tf.groupby('edge_gg').agg({
    'mean_importance': 'max',
    'median_importance': 'max'
}).reset_index()
net1_mean_imp = dict(zip(net1_agg['edge_gg'], net1_agg['mean_importance']))
net1_median_imp = dict(zip(net1_agg['edge_gg'], net1_agg['median_importance']))

# Net2: aggregate to gene-gene level (max across transcripts)
net2_agg = net2_tf.loc[net2_tf.groupby('edge_gg')['mean_importance'].idxmax()]
net2_mean_imp = dict(zip(net2_agg['edge_gg'], net2_agg['mean_importance']))
net2_median_imp = dict(zip(net2_agg['edge_gg'], net2_agg['median_importance']))
net2_best_tx = dict(zip(net2_agg['edge_gg'], net2_agg['source_transcript']))

all_edges_t1 = set(net1_mean_imp.keys()) | set(net2_mean_imp.keys())
# ---

# Cell 21
def categorize_source_resolution(edge_gg, net1_mean, net2_mean, net1_median, net2_median, net2_best_tx,
                                  r_iso=R_ISO, r_gene=R_GENE, r_eq=R_EQ, eps=EPSILON):
    """
    Categorize edge based on source resolution.
    Uses MEAN importance for ratio calculation (AlterNet 1.0 style).
    """
    S1_mean = net1_mean.get(edge_gg, 0)
    S2_mean = net2_mean.get(edge_gg, 0)
    S1_median = net1_median.get(edge_gg, 0)
    S2_median = net2_median.get(edge_gg, 0)
    best_tx = net2_best_tx.get(edge_gg, None)
    
    E1, E2 = S1_mean > 0, S2_mean > 0
    reg_gene, target_gene = edge_gg.split('|')
    
    # Use MEAN importance for fold-change calculation
    if E2 and not E1:
        category = 'source_isoform_specific'
    elif E1 and not E2:
        category = 'source_gene_specific'
    elif E1 and E2:
        ratio = S2_mean / (S1_mean + eps)
        if ratio >= r_iso:
            category = 'source_isoform_specific'
        elif ratio <= 1/r_gene:
            category = 'source_gene_specific'
        elif 1/r_eq <= ratio <= r_eq:
            category = 'source_equivalent'
        else:
            category = 'source_ambiguous'
    else:
        category = 'source_ambiguous'
    
    return {
        'regulator_gene': reg_gene, 'target_gene': target_gene, 'best_tx': best_tx,
        'S1_mean': S1_mean, 'S2_mean': S2_mean,
        'S1_median': S1_median, 'S2_median': S2_median,
        'E1': E1, 'E2': E2,
        'ratio': S2_mean / (S1_mean + eps) if E1 else np.inf,
        'source_category': category
    }
# ---

# Cell 22
set_a_rows = [categorize_source_resolution(e, net1_mean_imp, net2_mean_imp, 
                                             net1_median_imp, net2_median_imp, net2_best_tx) 
               for e in all_edges_t1]
set_a = pd.DataFrame(set_a_rows)

# Sort by median importance (AlterNet 1.0 style)
set_a['max_median'] = set_a[['S1_median', 'S2_median']].max(axis=1)
set_a = set_a.sort_values('max_median', ascending=False)

print(f"\nSet A: {len(set_a):,} rows")
print("\nCategory distribution:")
for cat, count in set_a['source_category'].value_counts().items():
    print(f"  {cat:30} {count:>8,} ({count/len(set_a)*100:>5.1f}%)")
# ---

# Cell 24
print(f"\nFiltered TF+SF edges from Net3: {len(net3_tfsf):,}")
# ---

# Cell 25
def compute_tfsf_evidence(reg_tx, target_tx, expr_df_indexed, usage_df_indexed, 
                          reliability_df_indexed, sf_expr_matrix, sample_cols):
    result = {
        'tf_rho': np.nan, 'tf_pval': np.nan,
        'tf_rho_conditional': np.nan, 'tf_pval_conditional': np.nan,
        'sf_rho': np.nan, 'sf_pval': np.nan, 
        'delta_usage': np.nan,
        'qc_ok': False, 'n_samples': 0
    }
    
    if reg_tx not in expr_df_indexed.index or target_tx not in expr_df_indexed.index:
        return result
    if target_tx not in usage_df_indexed.index:
        return result
    
    reg_expr = np.array(expr_df_indexed.loc[reg_tx].values, dtype=np.float64)
        target_expr = np.array(expr_df_indexed.loc[target_tx].values, dtype=np.float64)
        target_usage = np.array(usage_df_indexed.loc[target_tx].values, dtype=np.float64)
        target_reliable = np.array(reliability_df_indexed.loc[target_tx].values, dtype=np.float64)
        
        if reg_expr.ndim > 1: reg_expr = reg_expr[0]
        if target_expr.ndim > 1: target_expr = target_expr[0]
        if target_usage.ndim > 1: target_usage = target_usage[0]
        if target_reliable.ndim > 1: target_reliable = target_reliable[0]
        
        valid = (target_reliable > 0) & np.isfinite(reg_expr) & np.isfinite(target_expr) & np.isfinite(target_usage)
        if valid.sum() < 20:
            return result
        
        result['n_samples'] = int(valid.sum())
        result['qc_ok'] = True
        
        tf_rho, tf_pval = spearmanr(reg_expr[valid], target_expr[valid])
        result['tf_rho'] = tf_rho
        result['tf_pval'] = tf_pval
        
        if sf_expr_matrix is not None and len(sf_expr_matrix.columns) > 0:
            from sklearn.linear_model import LinearRegression
            valid_indices = np.where(valid)[0]
            sf_data = sf_expr_matrix.iloc[valid_indices].values
            sf_var = np.var(sf_data, axis=0)
            sf_cols_valid = sf_var > 1e-6
            if sf_cols_valid.sum() > 0:
                sf_data_filtered = sf_data[:, sf_cols_valid]
                lr = LinearRegression()
                lr.fit(sf_data_filtered, target_expr[valid])
                target_residuals = target_expr[valid] - lr.predict(sf_data_filtered)
                tf_rho_cond, tf_pval_cond = spearmanr(reg_expr[valid], target_residuals)
                result['tf_rho_conditional'] = tf_rho_cond
                result['tf_pval_conditional'] = tf_pval_cond
            else:
                result['tf_rho_conditional'] = tf_rho
                result['tf_pval_conditional'] = tf_pval
        else:
            result['tf_rho_conditional'] = tf_rho
            result['tf_pval_conditional'] = tf_pval
        
        sf_rho, sf_pval = spearmanr(reg_expr[valid], target_usage[valid])
        result['sf_rho'] = sf_rho
        result['sf_pval'] = sf_pval
        
        q25, q75 = np.percentile(reg_expr[valid], [25, 75])
        low_mask, high_mask = reg_expr[valid] <= q25, reg_expr[valid] >= q75
        if low_mask.sum() > 0 and high_mask.sum() > 0:
            result['delta_usage'] = np.mean(target_usage[valid][high_mask]) - np.mean(target_usage[valid][low_mask])
        else:
            result['delta_usage'] = 0
    return result
# ---

# Cell 26
if len(net3_tfsf) > 0:
    sf_tx_list = regulator_list[regulator_list['Regulator_type'] == 'SF']['Transcript stable ID'].tolist()
    sf_tx_in_expr = [tx for tx in sf_tx_list if tx in expr_df_indexed.index]
    if len(sf_tx_in_expr) > 0:
        sf_expr_matrix = expr_df_indexed.loc[sf_tx_in_expr].T
        sf_expr_matrix = sf_expr_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
    else:
        sf_expr_matrix = None
    
    tfsf_valid = net3_tfsf[
        net3_tfsf['source_transcript'].isin(expr_df_indexed.index) &
        net3_tfsf['target_transcript'].isin(expr_df_indexed.index)
    ].copy()
# ---

# Cell 27
if len(tfsf_valid) > 0:
    
    set_d_rows = []
    for idx, (_, row) in enumerate(tfsf_valid.iterrows()):
        if idx % 500 == 0:
        evidence = compute_tfsf_evidence(
            row['source_transcript'], row['target_transcript'],
            expr_df_indexed, usage_df_indexed, reliability_df_indexed,
            sf_expr_matrix, sample_cols
        )
        set_d_rows.append({
            'reg_tx': row['source_transcript'],
            'reg_gene': row['source_gene'],
            'target_tx': row['target_transcript'],
            'target_gene': row['target_gene'],
            'mean_importance': row['mean_importance'],
            'median_importance': row['median_importance'],
            **evidence
        })
    
    set_d_full = pd.DataFrame(set_d_rows)
    
    # FDR correction
    valid_tf = set_d_full['tf_pval_conditional'].notna() & (set_d_full['tf_pval_conditional'] > 0)
    valid_sf = set_d_full['sf_pval'].notna() & (set_d_full['sf_pval'] > 0)
    set_d_full['tf_q'] = np.nan
    set_d_full['sf_q'] = np.nan
    if valid_tf.sum() > 0:
        _, tf_qvals, _, _ = multipletests(set_d_full.loc[valid_tf, 'tf_pval_conditional'].values, method='fdr_bh')
        set_d_full.loc[valid_tf, 'tf_q'] = tf_qvals
    if valid_sf.sum() > 0:
        _, sf_qvals, _, _ = multipletests(set_d_full.loc[valid_sf, 'sf_pval'].values, method='fdr_bh')
        set_d_full.loc[valid_sf, 'sf_q'] = sf_qvals
    
    set_d_full['tf_evidence'] = ((set_d_full['tf_q'] <= Q_MIN) & (set_d_full['tf_rho_conditional'].abs() >= RHO_MIN)).fillna(False)
    set_d_full['sf_evidence'] = ((set_d_full['sf_q'] <= Q_MIN) & (set_d_full['sf_rho'].abs() >= RHO_MIN) & (set_d_full['delta_usage'].abs() >= DU_MIN)).fillna(False)
    
    def categorize_tfsf(row):
        if not row['qc_ok']:
            return 'tfsf_ambiguous'
        tf, sf = row['tf_evidence'], row['sf_evidence']
        if sf and not tf:
            return 'tfsf_sf_like'
        elif tf and not sf:
            return 'tfsf_tf_like'
        elif tf and sf:
            return 'tfsf_joint'
        else:
            return 'tfsf_ambiguous'
    
    set_d_full['tfsf_category'] = set_d_full.apply(categorize_tfsf, axis=1)
    
    # Sort by median importance
    set_d_full = set_d_full.sort_values('median_importance', ascending=False)
    
else:
    set_d_full = net3_tfsf[['source_transcript', 'source_gene', 'target_transcript', 'target_gene', 
                             'mean_importance', 'median_importance']].copy()
    set_d_full.columns = ['reg_tx', 'reg_gene', 'target_tx', 'target_gene', 'mean_importance', 'median_importance']
    set_d_full['tfsf_category'] = 'tfsf_ambiguous'
# ---

# Cell 28
print("\nSet D - Initial TF+SF Categories:")
for cat, count in set_d_full['tfsf_category'].value_counts().items():
    print(f"  {cat:25} {count:>8,} ({count/len(set_d_full)*100:>5.1f}%)")

# Split for downstream processing
tfsf_tf_like = set_d_full[set_d_full['tfsf_category'] == 'tfsf_tf_like'].copy()
tfsf_sf_like = set_d_full[set_d_full['tfsf_category'] == 'tfsf_sf_like'].copy()
tfsf_joint = set_d_full[set_d_full['tfsf_category'] == 'tfsf_joint'].copy()
tfsf_ambiguous = set_d_full[set_d_full['tfsf_category'] == 'tfsf_ambiguous'].copy()

print(f"\n  tfsf_tf_like: {len(tfsf_tf_like):,} (will be added to Set B)")
print(f"  tfsf_sf_like: {len(tfsf_sf_like):,} (will be added to Set C)")
print(f"  tfsf_joint: {len(tfsf_joint):,} (stays in Set D)")
print(f"  tfsf_ambiguous: {len(tfsf_ambiguous):,} (stays in Set D)")
# ---

# Cell 30
# Net2: ALL TF-like edges (TF + TF_SF) - they all act as TFs in gene-level target network
net2_for_t2 = net2_filtered[net2_filtered['reg_type'].isin(['TF', 'TF_SF'])].copy()

# Net3: TF only (TF_SF is handled via Set D)
net3_tf_only = net3_filtered[net3_filtered['reg_type'] == 'TF'].copy()

print(f"Net2 (TF + TF_SF acting as TF): {len(net2_for_t2):,}")
print(f"Net3 (TF only): {len(net3_tf_only):,}")
# ---

# Cell 31
# Prepare Net2: transcript-gene level
net2_for_t2['edge_tg'] = net2_for_t2['source_transcript'] + '|' + net2_for_t2['target_gene']
net2_tg_mean_imp = dict(zip(net2_for_t2['edge_tg'], net2_for_t2['mean_importance']))
net2_tg_median_imp = dict(zip(net2_for_t2['edge_tg'], net2_for_t2['median_importance']))

# Prepare Net3 TF-only: aggregate to transcript-gene level
net3_tf_only['edge_tg'] = net3_tf_only['source_transcript'] + '|' + net3_tf_only['target_gene']
net3_tg_agg = net3_tf_only.groupby('edge_tg').agg({
    'mean_importance': ['sum', 'max', 'count'],
    'median_importance': ['sum', 'max'],
    'source_gene': 'first',
    'target_gene': 'first',
    'target_transcript': lambda x: ','.join(sorted(set(x)))
}).reset_index()
net3_tg_agg.columns = ['edge_tg', 'S3_mean_sum', 'S3_mean_max', 'n_target_tx', 
                       'S3_median_sum', 'S3_median_max', 'source_gene', 'target_gene',
                       'target_tx_set']
net3_tg_agg['dominance'] = net3_tg_agg['S3_mean_max'] / (net3_tg_agg['S3_mean_sum'] + EPSILON)

net3_tg_mean_imp = dict(zip(net3_tg_agg['edge_tg'], net3_tg_agg['S3_mean_sum']))
net3_tg_mean_max = dict(zip(net3_tg_agg['edge_tg'], net3_tg_agg['S3_mean_max']))
net3_tg_median_imp = dict(zip(net3_tg_agg['edge_tg'], net3_tg_agg['S3_median_sum']))
net3_tg_dom = dict(zip(net3_tg_agg['edge_tg'], net3_tg_agg['dominance']))

# All edges for TF part of Set B
all_edges_t2_tf = set(net2_tg_mean_imp.keys()) | set(net3_tg_mean_imp.keys())
print(f"\nTF edges (Net2 TF+TF_SF vs Net3 TF): {len(all_edges_t2_tf):,}")
# ---

# Cell 32
def categorize_target_resolution(edge_tg, net2_mean, net3_mean, net3_max, net3_dom,
                                  net2_median, net3_median,
                                  r_iso=R_ISO, r_gene=R_GENE, r_eq=R_EQ, eps=EPSILON):
    """
    Categorize edge based on target resolution.
    Uses MEAN importance for ratio calculation (AlterNet 1.0 style).
    """
    S2_mean = net2_mean.get(edge_tg, 0)
    S3_mean_sum = net3_mean.get(edge_tg, 0)
    S3_mean_max = net3_max.get(edge_tg, 0)
    S2_median = net2_median.get(edge_tg, 0)
    S3_median = net3_median.get(edge_tg, 0)
    dominance = net3_dom.get(edge_tg, 0)
    
    E2, E3 = S2_mean > 0, S3_mean_sum > 0
    reg_tx, target_gene = edge_tg.split('|')
    reg_gene = tx2gene.get(reg_tx, None)
    
    # Use MEAN importance for fold-change calculation
    if E3 and not E2:
        category = 'target_isoform_specific'
    elif E2 and not E3:
        category = 'target_gene_specific'
    elif E2 and E3:
        ratio = S3_mean_sum / (S2_mean + eps)
        if ratio >= r_iso:
            category = 'target_isoform_specific'
        elif ratio <= 1/r_gene:
            category = 'target_gene_specific'
        elif 1/r_eq <= ratio <= r_eq:
            category = 'target_equivalent'
        else:
            category = 'target_ambiguous'
    else:
        category = 'target_ambiguous'
    
    return {
        'regulator_tx': reg_tx, 'regulator_gene': reg_gene, 'target_gene': target_gene,
        'S2_mean': S2_mean, 'S3_mean_sum': S3_mean_sum, 'S3_mean_max': S3_mean_max,
        'S2_median': S2_median, 'S3_median': S3_median,
        'dominance': dominance,
        'E2': E2, 'E3': E3,
        'ratio': S3_mean_sum / (S2_mean + eps) if E2 else np.inf,
        'target_category': category
    }
# ---

# Cell 33
# Categorize TF edges (from Net2 TF+TF_SF vs Net3 TF)
set_b_tf_rows = [categorize_target_resolution(e, net2_tg_mean_imp, net3_tg_mean_imp, 
                                                net3_tg_mean_max, net3_tg_dom,
                                                net2_tg_median_imp, net3_tg_median_imp) 
                  for e in all_edges_t2_tf]
set_b_tf = pd.DataFrame(set_b_tf_rows)
set_b_tf['reg_type'] = 'TF'

# Attach target_tx_set from Net3 aggregation
_net3_tx_map = dict(zip(net3_tg_agg['edge_tg'], net3_tg_agg['target_tx_set']))
set_b_tf['target_tx_set'] = set_b_tf.apply(
    lambda r: _net3_tx_map.get(r['regulator_tx'] + '|' + r['target_gene'], ''), axis=1)

print(f"Set B base (TF): {len(set_b_tf):,} rows")
# ---

# Cell 34
# Now add tfsf_tf_like edges
if len(tfsf_tf_like) > 0:
    
    # Get Net2 TF_SF edges (for comparison)
    net2_tfsf = net2_filtered[net2_filtered['reg_type'] == 'TF_SF'].copy()
    net2_tfsf['edge_tg'] = net2_tfsf['source_transcript'] + '|' + net2_tfsf['target_gene']
    net2_tfsf_mean_imp = dict(zip(net2_tfsf['edge_tg'], net2_tfsf['mean_importance']))
    net2_tfsf_median_imp = dict(zip(net2_tfsf['edge_tg'], net2_tfsf['median_importance']))
    
    # Get Net3 TF_SF edges aggregated (for the tf_like subset)
    net3_tfsf_agg = net3_tfsf.groupby(['source_transcript', 'target_gene']).agg({
        'mean_importance': ['sum', 'max', 'count'],
        'median_importance': ['sum', 'max'],
        'source_gene': 'first',
        'target_transcript': lambda x: ','.join(sorted(set(x)))
    }).reset_index()
    net3_tfsf_agg.columns = ['source_transcript', 'target_gene', 'S3_mean_sum', 'S3_mean_max', 'n_target_tx',
                            'S3_median_sum', 'S3_median_max', 'source_gene', 'target_tx_set']
    net3_tfsf_agg['edge_tg'] = net3_tfsf_agg['source_transcript'] + '|' + net3_tfsf_agg['target_gene']
    net3_tfsf_agg['dominance'] = net3_tfsf_agg['S3_mean_max'] / (net3_tfsf_agg['S3_mean_sum'] + EPSILON)
    
    tfsf_tg_mean_imp = dict(zip(net3_tfsf_agg['edge_tg'], net3_tfsf_agg['S3_mean_sum']))
    tfsf_tg_mean_max = dict(zip(net3_tfsf_agg['edge_tg'], net3_tfsf_agg['S3_mean_max']))
    tfsf_tg_median_imp = dict(zip(net3_tfsf_agg['edge_tg'], net3_tfsf_agg['S3_median_sum']))
    tfsf_tg_dom = dict(zip(net3_tfsf_agg['edge_tg'], net3_tfsf_agg['dominance']))
    
    # Get unique edges from tfsf_tf_like
    tfsf_tf_like['edge_tg'] = tfsf_tf_like['reg_tx'] + '|' + tfsf_tf_like['target_gene']
    tfsf_tf_like_edges = set(tfsf_tf_like['edge_tg'].unique())
    
    # Categorize tfsf_tf_like edges
    set_b_tfsf_rows = [categorize_target_resolution(e, net2_tfsf_mean_imp, tfsf_tg_mean_imp, 
                                                      tfsf_tg_mean_max, tfsf_tg_dom,
                                                      net2_tfsf_median_imp, tfsf_tg_median_imp) 
                        for e in tfsf_tf_like_edges]
    set_b_tfsf = pd.DataFrame(set_b_tfsf_rows)
    set_b_tfsf['reg_type'] = 'TF_SF'
    
    # Attach target_tx_set
    _tfsf_tx_map = dict(zip(
        net3_tfsf_agg['source_transcript'] + '|' + net3_tfsf_agg['target_gene'],
        net3_tfsf_agg['target_tx_set']))
    set_b_tfsf['target_tx_set'] = set_b_tfsf.apply(
        lambda r: _tfsf_tx_map.get(r['regulator_tx'] + '|' + r['target_gene'], ''), axis=1)
    
else:
    set_b_tfsf = pd.DataFrame()
# ---

# Cell 35
# Combine into final Set B
set_b = pd.concat([set_b_tf, set_b_tfsf], ignore_index=True)

# Sort by median importance (AlterNet 1.0 style)
set_b['max_median'] = set_b[['S2_median', 'S3_median']].max(axis=1)
set_b = set_b.sort_values('max_median', ascending=False)

print(f"\nFinal Set B: {len(set_b):,} rows")
print(f"  - TF edges (from Net2 TF+TF_SF vs Net3 TF): {len(set_b_tf):,}")
print(f"  - TF_SF (tf_like from Set D): {len(set_b_tfsf) if len(set_b_tfsf) > 0 else 0:,}")
print("\nCategory distribution:")
for cat, count in set_b['target_category'].value_counts().items():
    print(f"  {cat:30} {count:>8,} ({count/len(set_b)*100:>5.1f}%)")
print("\nBy regulator type:")
print(set_b.groupby(['reg_type', 'target_category']).size().unstack(fill_value=0))
# ---

# Cell 37
# Set B Unpacked: transcript-resolved edges from Net3

# Build unpacked table from Net3 individual edges

# Source 1: TF edges from net3_tf_only
net3_tf_unpack = net3_tf_only[['source_transcript', 'source_gene', 'target_transcript',
                                'target_gene', 'mean_importance', 'median_importance', 
                                'frequency']].copy()
net3_tf_unpack.columns = ['regulator_tx', 'regulator_gene', 'target_tx', 'target_gene',
                           'net3_mean_importance', 'net3_median_importance', 'net3_frequency']
net3_tf_unpack['reg_type'] = 'TF'

# Source 2: tfsf_tf_like edges from net3_tfsf
if len(tfsf_tf_like) > 0:
    tfsf_tf_like_regtx_set = set(tfsf_tf_like['reg_tx'].unique())
    net3_tfsf_unpack = net3_tfsf[
        net3_tfsf['source_transcript'].isin(tfsf_tf_like_regtx_set)
    ][['source_transcript', 'source_gene', 'target_transcript',
       'target_gene', 'mean_importance', 'median_importance', 
       'frequency']].copy()
    net3_tfsf_unpack.columns = ['regulator_tx', 'regulator_gene', 'target_tx', 'target_gene',
                                 'net3_mean_importance', 'net3_median_importance', 'net3_frequency']
    net3_tfsf_unpack['reg_type'] = 'TF_SF'
    set_b_unpacked = pd.concat([net3_tf_unpack, net3_tfsf_unpack], ignore_index=True)
else:
    set_b_unpacked = net3_tf_unpack.copy()

# Join with Set B categorization
# Create join key
set_b_unpacked['edge_tg'] = set_b_unpacked['regulator_tx'] + '|' + set_b_unpacked['target_gene']
set_b_key = set_b[['regulator_tx', 'target_gene', 'target_category', 'dominance',
                      'S3_mean_sum', 'target_tx_set']].copy()
set_b_key['edge_tg'] = set_b_key['regulator_tx'] + '|' + set_b_key['target_gene']

# Merge: only keep edges whose parent exists in Set B
set_b_unpacked = set_b_unpacked.merge(
    set_b_key[['edge_tg', 'target_category']],
    on='edge_tg', how='inner'
)

# Add target_rank_within_edge
set_b_unpacked['target_rank_within_edge'] = set_b_unpacked.groupby('edge_tg')[
    'net3_mean_importance'].rank(ascending=False, method='min').astype(int)

# Sort by importance
set_b_unpacked = set_b_unpacked.sort_values('net3_median_importance', ascending=False)

# Drop join key
set_b_unpacked = set_b_unpacked.drop(columns=['edge_tg'])

print(f"Set B Unpacked: {len(set_b_unpacked):,} transcript-level edges")
print(f"  From TF edges: {(set_b_unpacked['reg_type'] == 'TF').sum():,}")
print(f"  From TF_SF (tf_like): {(set_b_unpacked['reg_type'] == 'TF_SF').sum():,}")
print(f"\nCategory distribution (inherited):")
for cat, count in set_b_unpacked['target_category'].value_counts().items():
    print(f"  {cat:30} {count:>8,} ({count/len(set_b_unpacked)*100:>5.1f}%)")
print(f"\nRank distribution:")
print(f"  Rank 1 (dominant target tx): {(set_b_unpacked['target_rank_within_edge'] == 1).sum():,}")
print(f"  Rank 2+: {(set_b_unpacked['target_rank_within_edge'] > 1).sum():,}")

# ---

# Cell 38
# Add presentation columns to Set B

def get_presentation_cols_t2(row):
    """Derive target_tx_resolved, target_tx_dominant, n_target_tx from target_tx_set."""
    tx_set_str = row.get('target_tx_set', '')
    if pd.isna(tx_set_str) or tx_set_str == '':
        return pd.Series({'n_target_tx': 0, 'target_tx_resolved': '', 'target_tx_dominant': ''})
    
    tx_list = [t.strip() for t in str(tx_set_str).split(',') if t.strip()]
    n = len(tx_list)
    
    # Resolved: if exactly 1 transcript
    resolved = tx_list[0] if n == 1 else ''
    
    return pd.Series({'n_target_tx': n, 'target_tx_resolved': resolved, 'target_tx_dominant': ''})

# Apply basic columns
pres_cols = set_b.apply(get_presentation_cols_t2, axis=1)
set_b['n_target_tx'] = pres_cols['n_target_tx'].astype(int)
set_b['target_tx_resolved'] = pres_cols['target_tx_resolved']

# For target_tx_dominant: find the argmax target_tx by Net3 importance
# Use the unpacked table which has per-edge per-target importance
_dominant_map = set_b_unpacked[
    set_b_unpacked['target_rank_within_edge'] == 1
].set_index(
    set_b_unpacked[set_b_unpacked['target_rank_within_edge'] == 1].apply(
        lambda r: r['regulator_tx'] + '|' + r['target_gene'], axis=1)
)['target_tx'].to_dict()

set_b['target_tx_dominant'] = set_b.apply(
    lambda r: _dominant_map.get(r['regulator_tx'] + '|' + r['target_gene'], ''), axis=1)

print(f"  n_target_tx == 0 (gene_specific, no Net3 edges): {(set_b['n_target_tx'] == 0).sum():,}")
print(f"  n_target_tx == 1 (resolved): {(set_b['n_target_tx'] == 1).sum():,}")
print(f"  n_target_tx >= 2 (multi): {(set_b['n_target_tx'] >= 2).sum():,}")
print(f"  target_tx_dominant non-empty: {(set_b['target_tx_dominant'] != '').sum():,}")

# ---

# Cell 40
print(f"\nFiltered SF edges from Net3: {len(net3_sf):,}")
# ---

# Cell 41
# Aggregate SF edges to (SF_tx, target_gene) level
sf_edges = net3_sf.groupby(['source_transcript', 'source_gene', 'target_gene']).agg({
    'mean_importance': ['sum', 'max', 'count'],
    'median_importance': ['sum', 'max'],
    'target_transcript': lambda x: set(x)
}).reset_index()
sf_edges.columns = ['sf_tx', 'sf_gene', 'target_gene', 'mean_importance_sum', 'mean_importance_max', 
                    'n_target_tx', 'median_importance_sum', 'median_importance_max', 'target_tx_set']

print(f"Aggregated SF edges (sf_tx, target_gene): {len(sf_edges):,}")
# ---

# Cell 42
def compute_sf_splicing_evidence(sf_tx, target_gene, target_tx_set, mean_importance_sum, mean_importance_max,
                                  median_importance_sum, median_importance_max,
                                  sf_gene, n_target_tx, sf_expr_dict):
    result = {
        'sf_tx': sf_tx, 'sf_gene': sf_gene, 'target_gene': target_gene,
        'target_tx_set': ','.join(sorted(target_tx_set)) if isinstance(target_tx_set, set) else str(target_tx_set),
        'mean_importance_sum': mean_importance_sum, 'mean_importance_max': mean_importance_max,
        'median_importance_sum': median_importance_sum, 'median_importance_max': median_importance_max,
        'n_target_tx_in_net': n_target_tx,
        'n_isoforms_in_gene': 0, 'n_transcripts_tested': 0,
        'n_sig': 0, 'push_pull': False, 'dominance': 0.0,
        'usage_reliable': False, 'best_rho': 0.0, 'best_tx': None,
        'sf_category': 'sf_ambiguous'
    }
    
    if sf_tx not in sf_expr_dict:
        return result
    
    sf_vals = sf_expr_dict[sf_tx]
    all_gene_txs = gene_to_all_transcripts.get(target_gene, [])
    all_gene_txs_in_data = [tx for tx in all_gene_txs if tx in usage_df_indexed.index]
    result['n_isoforms_in_gene'] = len(all_gene_txs_in_data)
    
    if len(all_gene_txs_in_data) < MIN_ISOFORMS_FOR_SPLICING:
        return result
    
    txs_to_test = set(tx for tx in target_tx_set if tx in usage_df_indexed.index)
    gene_txs_sorted = sorted(all_gene_txs_in_data, key=lambda tx: transcript_mean_expr.get(tx, 0), reverse=True)
    for tx in gene_txs_sorted[:TOP_M_EXPRESSED]:
        txs_to_test.add(tx)
    txs_to_test = list(txs_to_test)
    
    if len(txs_to_test) == 0:
        return result
    
    result['n_transcripts_tested'] = len(txs_to_test)
    result['usage_reliable'] = True
    
    correlations = []
    for tx in txs_to_test:
            tx_usage = np.array(usage_df_indexed.loc[tx].values, dtype=np.float64)
            tx_reliable = np.array(reliability_df_indexed.loc[tx].values, dtype=np.float64)
            if tx_usage.ndim > 1: tx_usage = tx_usage[0]
            if tx_reliable.ndim > 1: tx_reliable = tx_reliable[0]
            
            valid = (tx_reliable > 0) & np.isfinite(sf_vals) & np.isfinite(tx_usage)
            if valid.sum() < 10: continue
            
            rho, pval = spearmanr(sf_vals[valid], tx_usage[valid])
            if not np.isfinite(rho): continue
            
            q25, q75 = np.percentile(sf_vals[valid], [25, 75])
            low_mask, high_mask = sf_vals[valid] <= q25, sf_vals[valid] >= q75
            delta_usage = 0
            if low_mask.sum() > 0 and high_mask.sum() > 0:
                delta_usage = np.mean(tx_usage[valid][high_mask]) - np.mean(tx_usage[valid][low_mask])
            
            correlations.append({'tx': tx, 'rho': rho, 'pval': pval, 'delta_usage': delta_usage})
    
    if len(correlations) == 0:
        result['sf_category'] = 'sf_expression_associated'
        return result
    
    _, qvals, _, _ = multipletests([c['pval'] for c in correlations], method='fdr_bh')
    for i, c in enumerate(correlations): c['qval'] = qvals[i]
    
    sig_corrs = [c for c in correlations if c['qval'] <= Q_MIN and abs(c['rho']) >= RHO_MIN and abs(c['delta_usage']) >= DU_MIN]
    n_sig = len(sig_corrs)
    result['n_sig'] = n_sig
    
    if n_sig >= 2:
        rhos = [c['rho'] for c in sig_corrs]
        result['push_pull'] = (max(rhos) > 0 and min(rhos) < 0)
    
    if n_sig > 0:
        abs_rhos = [abs(c['rho']) for c in sig_corrs]
        result['dominance'] = max(abs_rhos) / (sum(abs_rhos) + EPSILON)
        best_idx = np.argmax(abs_rhos)
        result['best_rho'] = sig_corrs[best_idx]['rho']
        result['best_tx'] = sig_corrs[best_idx]['tx']
    
    if n_sig >= 2:
        if result['push_pull'] or result['dominance'] >= DOM_MIN:  # DOM_MIN = 0.5
            result['sf_category'] = 'sf_splicing_supported_specific'
        else:
            result['sf_category'] = 'sf_splicing_supported_broad'
    else:
        result['sf_category'] = 'sf_expression_associated'
    
    return result
# ---

# Cell 43
# Build SF expression dict
sf_transcripts = sf_edges['sf_tx'].unique()
sf_expr = transcript_tpm[transcript_tpm['transcript_id'].isin(sf_transcripts)].copy()
sf_expr = sf_expr.set_index('transcript_id')[sample_cols]
sf_expr_dict = {tx: np.array(sf_expr.loc[tx].values, dtype=np.float64) for tx in sf_expr.index}
# ---

# Cell 44
if len(sf_expr_dict) > 0:
    
    set_c_sf_rows = []
    for idx, row in sf_edges.iterrows():
        if idx % 1000 == 0:
        result = compute_sf_splicing_evidence(
            row['sf_tx'], row['target_gene'], row['target_tx_set'],
            row['mean_importance_sum'], row['mean_importance_max'],
            row['median_importance_sum'], row['median_importance_max'],
            row['sf_gene'], row['n_target_tx'], sf_expr_dict
        )
        set_c_sf_rows.append(result)
    
    set_c_sf = pd.DataFrame(set_c_sf_rows)
    set_c_sf['reg_type'] = 'SF'
    print(f"\nSet C SF-only: {len(set_c_sf):,} rows")
else:
    set_c_sf = sf_edges.copy()
    set_c_sf['sf_category'] = 'sf_expression_associated'
    set_c_sf['reg_type'] = 'SF'
# ---

# Cell 45
# Now add tfsf_sf_like edges
if len(tfsf_sf_like) > 0:
    
    # Aggregate tfsf_sf_like to (reg_tx, target_gene) level
    tfsf_sf_agg = net3_tfsf.groupby(['source_transcript', 'source_gene', 'target_gene']).agg({
        'mean_importance': ['sum', 'max', 'count'],
        'median_importance': ['sum', 'max'],
        'target_transcript': lambda x: set(x)
    }).reset_index()
    tfsf_sf_agg.columns = ['sf_tx', 'sf_gene', 'target_gene', 'mean_importance_sum', 'mean_importance_max',
                          'n_target_tx', 'median_importance_sum', 'median_importance_max', 'target_tx_set']
    
    # Filter to sf_like edges only
    tfsf_sf_like_regtx = set(tfsf_sf_like['reg_tx'].unique())
    tfsf_sf_agg_filtered = tfsf_sf_agg[tfsf_sf_agg['sf_tx'].isin(tfsf_sf_like_regtx)].copy()
    
    # Add TF_SF transcripts to sf_expr_dict
    tfsf_transcripts = tfsf_sf_agg_filtered['sf_tx'].unique()
    tfsf_expr = transcript_tpm[transcript_tpm['transcript_id'].isin(tfsf_transcripts)].copy()
    tfsf_expr = tfsf_expr.set_index('transcript_id')[sample_cols]
    for tx in tfsf_expr.index:
        sf_expr_dict[tx] = np.array(tfsf_expr.loc[tx].values, dtype=np.float64)
    
    # Categorize
    set_c_tfsf_rows = []
    for idx, row in tfsf_sf_agg_filtered.iterrows():
        result = compute_sf_splicing_evidence(
            row['sf_tx'], row['target_gene'], row['target_tx_set'],
            row['mean_importance_sum'], row['mean_importance_max'],
            row['median_importance_sum'], row['median_importance_max'],
            row['sf_gene'], row['n_target_tx'], sf_expr_dict
        )
        set_c_tfsf_rows.append(result)
    
    set_c_tfsf = pd.DataFrame(set_c_tfsf_rows)
    set_c_tfsf['reg_type'] = 'TF_SF'
else:
    set_c_tfsf = pd.DataFrame()
# ---

# Cell 46
# Combine into final Set C
set_c = pd.concat([set_c_sf, set_c_tfsf], ignore_index=True)

# Sort by median importance (AlterNet 1.0 style)
set_c = set_c.sort_values('median_importance_sum', ascending=False)

print(f"\nFinal Set C: {len(set_c):,} rows")
print(f"  - SF edges: {len(set_c_sf):,}")
print(f"  - TF_SF (sf_like from Set D): {len(set_c_tfsf) if len(set_c_tfsf) > 0 else 0:,}")
print("\nCategory distribution:")
for cat, count in set_c['sf_category'].value_counts().items():
    print(f"  {cat:35} {count:>8,} ({count/len(set_c)*100:>5.1f}%)")
print("\nBy regulator type:")
print(set_c.groupby(['reg_type', 'sf_category']).size().unstack(fill_value=0))
# ---

# Cell 48
# Set C Unpacked: transcript-resolved SF edges from Net3

# Build unpacked table from Net3 individual edges

# Source 1: SF edges from net3_sf
net3_sf_unpack = net3_sf[['source_transcript', 'source_gene', 'target_transcript',
                           'target_gene', 'mean_importance', 'median_importance',
                           'frequency']].copy()
net3_sf_unpack.columns = ['sf_tx', 'sf_gene', 'target_tx', 'target_gene',
                           'net3_mean_importance', 'net3_median_importance', 'net3_frequency']
net3_sf_unpack['reg_type'] = 'SF'

# Source 2: tfsf_sf_like edges from net3_tfsf
if len(tfsf_sf_like) > 0:
    tfsf_sf_like_regtx_set = set(tfsf_sf_like['reg_tx'].unique())
    net3_tfsf_sf_unpack = net3_tfsf[
        net3_tfsf['source_transcript'].isin(tfsf_sf_like_regtx_set)
    ][['source_transcript', 'source_gene', 'target_transcript',
       'target_gene', 'mean_importance', 'median_importance',
       'frequency']].copy()
    net3_tfsf_sf_unpack.columns = ['sf_tx', 'sf_gene', 'target_tx', 'target_gene',
                                    'net3_mean_importance', 'net3_median_importance', 'net3_frequency']
    net3_tfsf_sf_unpack['reg_type'] = 'TF_SF'
    set_c_unpacked = pd.concat([net3_sf_unpack, net3_tfsf_sf_unpack], ignore_index=True)
else:
    set_c_unpacked = net3_sf_unpack.copy()

# Join with Set C categorization
set_c_unpacked['edge_sg'] = set_c_unpacked['sf_tx'] + '|' + set_c_unpacked['target_gene']
set_c_key = set_c[['sf_tx', 'target_gene', 'sf_category']].copy()
set_c_key['edge_sg'] = set_c_key['sf_tx'] + '|' + set_c_key['target_gene']

set_c_unpacked = set_c_unpacked.merge(
    set_c_key[['edge_sg', 'sf_category']],
    on='edge_sg', how='inner'
)

# Add target_rank_within_edge
set_c_unpacked['target_rank_within_edge'] = set_c_unpacked.groupby('edge_sg')[
    'net3_mean_importance'].rank(ascending=False, method='min').astype(int)

# Sort
set_c_unpacked = set_c_unpacked.sort_values('net3_median_importance', ascending=False)
set_c_unpacked = set_c_unpacked.drop(columns=['edge_sg'])

print(f"Set C Unpacked: {len(set_c_unpacked):,} transcript-level edges")
print(f"  From SF edges: {(set_c_unpacked['reg_type'] == 'SF').sum():,}")
print(f"  From TF_SF (sf_like): {(set_c_unpacked['reg_type'] == 'TF_SF').sum():,}")
print(f"\nCategory distribution (inherited):")
for cat, count in set_c_unpacked['sf_category'].value_counts().items():
    print(f"  {cat:35} {count:>8,} ({count/len(set_c_unpacked)*100:>5.1f}%)")

# ---

# Cell 49
# Add presentation columns to Set C

def get_presentation_cols_t3(row):
    tx_set_str = row.get('target_tx_set', '')
    if pd.isna(tx_set_str) or tx_set_str == '':
        return pd.Series({'n_target_tx': 0, 'target_tx_resolved': '', 'target_tx_dominant': ''})
    
    # target_tx_set may be a string or already processed
    tx_list = [t.strip() for t in str(tx_set_str).split(',') if t.strip()]
    n = len(tx_list)
    resolved = tx_list[0] if n == 1 else ''
    
    return pd.Series({'n_target_tx': n, 'target_tx_resolved': resolved, 'target_tx_dominant': ''})

pres_cols3 = set_c.apply(get_presentation_cols_t3, axis=1)
set_c['n_target_tx'] = pres_cols3['n_target_tx'].astype(int)
set_c['target_tx_resolved'] = pres_cols3['target_tx_resolved']

# For target_tx_dominant: use best_tx from splicing correlation if available,
# otherwise use importance-based argmax from unpacked table
_dominant_map3 = set_c_unpacked[
    set_c_unpacked['target_rank_within_edge'] == 1
].copy()
_dominant_map3['edge_sg'] = _dominant_map3['sf_tx'] + '|' + _dominant_map3['target_gene']
_dominant_dict3 = dict(zip(_dominant_map3['edge_sg'], _dominant_map3['target_tx']))

def get_dominant_t3(row):
    # Prefer best_tx from splicing correlation (it has biological meaning)
    if pd.notna(row.get('best_tx', None)) and row.get('best_tx', '') != '':
        return row['best_tx']
    # Fallback to importance-based argmax
    key = row['sf_tx'] + '|' + row['target_gene']
    return _dominant_dict3.get(key, '')

set_c['target_tx_dominant'] = set_c.apply(get_dominant_t3, axis=1)

print(f"  n_target_tx == 1 (resolved): {(set_c['n_target_tx'] == 1).sum():,}")
print(f"  n_target_tx >= 2 (multi): {(set_c['n_target_tx'] >= 2).sum():,}")
print(f"  target_tx_dominant non-empty: {(set_c['target_tx_dominant'] != '').sum():,}")
print(f"  target_tx_dominant from best_tx: {set_c['best_tx'].notna().sum() - (set_c['best_tx'] == '').sum():,}")

# ---

# Cell 51
# Final Set D contains only joint and ambiguous
set_d = pd.concat([tfsf_joint, tfsf_ambiguous], ignore_index=True)

# Sort by median importance
set_d = set_d.sort_values('median_importance', ascending=False)

print(f"\nFinal Set D: {len(set_d):,} rows")
print("\nCategory distribution:")
for cat, count in set_d['tfsf_category'].value_counts().items():
    print(f"  {cat:25} {count:>8,} ({count/len(set_d)*100:>5.1f}%)")
# ---

# Cell 52
# Save tables

set_a.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_a_source_resolution.tsv"), sep='\t', index=False)
set_b.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_b_target_resolution.tsv"), sep='\t', index=False)
set_c.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_c_sf_splicing.tsv"), sep='\t', index=False)
set_d.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_d_tfsf_joint_ambiguous.tsv"), sep='\t', index=False)

# Save unpacked transcript-level tables
set_b_unpacked.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_b_unpacked_txlevel.tsv"), sep='\t', index=False)
set_c_unpacked.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_set_c_unpacked_txlevel.tsv"), sep='\t', index=False)

net1_filtered.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_net1_filtered.tsv"), sep='\t', index=False)
net2_filtered.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_net2_filtered.tsv"), sep='\t', index=False)
net3_filtered.to_csv(op.join(PART2_RESULTS, f"{CONDITION}_net3_filtered.tsv"), sep='\t', index=False)

# ---

# Cell 53
# Save summary YAML
summary = {
    'condition': CONDITION,
    'version': 'v2_alternet1_style',
    'importance_usage': {
        'filtering': 'median_importance',
        'fold_change_ratio': 'mean_importance',
        'sorting_ranking': 'median_importance'
    },
    'design_notes': {
        'set_b': 'Net2 (TF+TF_SF all act as TF) vs Net3 (TF only) + tfsf_tf_like from Set_D',
        'set_c': 'Net3 (SF only) + tfsf_sf_like from Set_D',
        'set_d': 'Only tfsf_joint and tfsf_ambiguous remain',
        'set_b_unpacked': 'Individual tx->tx edges from Net3, with target_category inherited from Set B',
        'set_c_unpacked': 'Individual sf_tx->tx edges from Net3, with sf_category inherited from Set C'
    },
    'filtered_networks': {
        'net1': int(len(net1_filtered)),
        'net2': int(len(net2_filtered)),
        'net3': int(len(net3_filtered)),
        'net3_tf': int(len(net3_tf)),
        'net3_sf': int(len(net3_sf)),
        'net3_tfsf': int(len(net3_tfsf))
    },
    'set_a': {k: int(v) for k, v in set_a['source_category'].value_counts().to_dict().items()},
    'set_b': {
        'total': int(len(set_b)),
        'tf_edges': int((set_b['reg_type'] == 'TF').sum()),
        'tfsf_tf_like_edges': int((set_b['reg_type'] == 'TF_SF').sum()),
        'categories': {k: int(v) for k, v in set_b['target_category'].value_counts().to_dict().items()}
    },
    'set_b_unpacked': {
        'total': int(len(set_b_unpacked)),
        'rank1_edges': int((set_b_unpacked['target_rank_within_edge'] == 1).sum())
    },
    'set_c': {
        'total': int(len(set_c)),
        'sf_edges': int((set_c['reg_type'] == 'SF').sum()),
        'tfsf_sf_like_edges': int((set_c['reg_type'] == 'TF_SF').sum()),
        'categories': {k: int(v) for k, v in set_c['sf_category'].value_counts().to_dict().items()}
    },
    'set_c_unpacked': {
        'total': int(len(set_c_unpacked)),
        'rank1_edges': int((set_c_unpacked['target_rank_within_edge'] == 1).sum())
    },
    'set_d': {
        'total': int(len(set_d)),
        'categories': {k: int(v) for k, v in set_d['tfsf_category'].value_counts().to_dict().items()}
    }
}

with open(op.join(PART2_RESULTS, f"{CONDITION}_part2_summary_v2.yaml"), 'w') as f:
    yaml.dump(summary, f, default_flow_style=False)

# ---