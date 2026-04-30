
import pandas as pd
import numpy as np
import os
import os.path as op
import yaml
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression


def categorize_source_resolution(edge_gg, net1_mean, net2_mean, net1_median, net2_median, net2_best_tx,
                                  r_iso=1.5, r_gene=1.5, r_eq=1.2, eps=1e-6):
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
        'edge_key': edge_gg,
        'regulator_gene': reg_gene, 'target_gene': target_gene, 'best_tx': best_tx,
        'S1_mean': S1_mean, 'S2_mean': S2_mean,
        'S1_median': S1_median, 'S2_median': S2_median,
        'E1': E1, 'E2': E2,
        'ratio': S2_mean / (S1_mean + eps) if E1 else np.inf,
        'source_category': category
    }

def canonical_vs_source_as(net1_tf,net2_tf):

    net1_tf['edge_gg'] = net1_tf['source_gene'] + '_' + net1_tf['target_gene']
    net2_tf['edge_gg'] = net2_tf['source_gene'] + '_' + net2_tf['target_gene']

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
    

    set_a_rows = [categorize_source_resolution(e, net1_mean_imp, net2_mean_imp, net1_median_imp, net2_median_imp, net2_best_tx) for e in all_edges_t1]
    set_a = pd.DataFrame(set_a_rows)

    # Sort by median importance (AlterNet 1.0 style)
    set_a['max_median'] = set_a[['S1_median', 'S2_median']].max(axis=1)
    set_a = set_a.sort_values('max_median', ascending=False)

    print(f"\nSet A: {len(set_a):,} rows")
    print("\nCategory distribution:")
    for cat, count in set_a['source_category'].value_counts().items():
        print(f"  {cat:30} {count:>8,} ({count/len(set_a)*100:>5.1f}%)")

    return set_a


def calculate_transcript_usage(df, gene_col='gene_id', transcript_col='transcript_id', 
                               min_gene_tpm=1.0, epsilon=1e-8):
    """
    Calculates transcript usage (fraction of gene expression) and reliability masks.
    
    Returns:
        usage_df: DataFrame of transcript usage (0.0 to 1.0)
        reliability_df: Binary mask (1.0 if gene expression >= min_gene_tpm)
        metadata: Dictionary containing gene-to-transcript mapping and mean expressions
    """
    # 1. Identify numeric sample columns
    id_cols = {gene_col, transcript_col}
    sample_cols = [c for c in df.columns 
                   if c not in id_cols and pd.api.types.is_numeric_dtype(df[c])]
    
    # 2. Extract numeric data and gene grouping
    tpm_numeric = df[sample_cols].astype(float)
    genes = df[gene_col].values
    
    # 3. Calculate Gene Totals (optimized groupby-transform)
    # Using a temporary series for grouping to avoid dataframe copies
    gene_totals = tpm_numeric.groupby(genes).transform('sum')
    
    # 4. Calculate Usage and Reliability
    usage_values = tpm_numeric / (gene_totals + epsilon)
    reliability_values = (gene_totals >= min_gene_tpm).astype(float)
    
    # 5. Construct output DataFrames
    usage_df = pd.concat([df[[transcript_col, gene_col]], usage_values], axis=1)
    
    # Reliability keeps transcript_id as index for easy alignment
    reliability_df = reliability_values.copy()
    reliability_df.index = df[transcript_col].values
    
    # 6. Metadata generation
    gene_to_transcripts = df.groupby(gene_col)[transcript_col].apply(list).to_dict()
    mean_expr = tpm_numeric.mean(axis=1).to_dict()
    # Mapping mean_expr to transcript IDs
    mean_expr_map = dict(zip(df[transcript_col], mean_expr))

    metadata = {
        "gene_map": gene_to_transcripts,
        "mean_expression": mean_expr_map,
        "sample_columns": sample_cols
    }
    
    return usage_df, reliability_df, metadata


def compute_tfsf_evidence(reg_tx, target_tx, expr_df_indexed, usage_df_indexed, 
                          reliability_df_indexed, sf_expr_matrix):
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
    
    #Compute correlation between the expression of source and target
    tf_rho, tf_pval = spearmanr(reg_expr[valid], target_expr[valid])
    result['tf_rho'] = tf_rho
    result['tf_pval'] = tf_pval
    
    if sf_expr_matrix is not None and len(sf_expr_matrix.columns) > 0:
        from sklearn.linear_model import LinearRegression
        valid_indices = np.where(valid)[0]
        sf_data = sf_expr_matrix.iloc[valid_indices].values
        sf_var = np.var(sf_data, axis=0)
        sf_cols_valid = sf_var > 1e-6
        # Regress target on the splicefactor data and get conditional correlation
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
    
    # Compute correlation between the target usage >> Splice evidence
    sf_rho, sf_pval = spearmanr(reg_expr[valid], target_usage[valid])
    result['sf_rho'] = sf_rho
    result['sf_pval'] = sf_pval
    
    ## Shift in expression
    q25, q75 = np.percentile(reg_expr[valid], [25, 75])
    low_mask, high_mask = reg_expr[valid] <= q25, reg_expr[valid] >= q75
    if low_mask.sum() > 0 and high_mask.sum() > 0:
        result['delta_usage'] = np.mean(target_usage[valid][high_mask]) - np.mean(target_usage[valid][low_mask])
    else:
        result['delta_usage'] = 0
    return result



def tf_sf_disambigouation_fully_as_aware(net3_tfsf, regulator_list, transcript_data, usage_df, reliability_df, sample_cols,
                                         RHO_MIN = 0.3,  Q_MIN = 0.05, DU_MIN = 0.1):

    # Step 1 get the data containing the sf and tf
    if len(net3_tfsf) > 0:
        sf_tx_list = regulator_list[regulator_list['Regulator_type'] == 'SF']['Transcript stable ID'].tolist()
        sf_tx_in_expr = [tx for tx in sf_tx_list if tx in transcript_data.index]
        if len(sf_tx_in_expr) > 0:
            sf_expr_matrix = transcript_data.loc[sf_tx_in_expr].T
            sf_expr_matrix = sf_expr_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
        else:
            sf_expr_matrix = None
        
        tfsf_valid = net3_tfsf[
            net3_tfsf['source_transcript'].isin(transcript_data.index) &
            net3_tfsf['target_transcript'].isin(transcript_data.index)
        ].copy()


    if len(tfsf_valid) > 0:
        
        set_d_rows = []
        for idx, (_, row) in enumerate(tfsf_valid.iterrows()):
            if idx % 500 == 0:
                evidence = compute_tfsf_evidence(
                    row['source_transcript'], row['target_transcript'],
                    transcript_data, usage_df, reliability_df,
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
        print('Warning falling back to classifying everything ambigouus')
        set_d_full = net3_tfsf[['source_transcript', 'source_gene', 'target_transcript', 'target_gene', 
                                'mean_importance', 'median_importance']].copy()
        set_d_full.columns = ['reg_tx', 'reg_gene', 'target_tx', 'target_gene', 'mean_importance', 'median_importance']
        set_d_full['tfsf_category'] = 'tfsf_ambiguous'
    
    return set_d_full
