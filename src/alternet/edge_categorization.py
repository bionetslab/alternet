
import pandas as pd
import numpy as np
import os
import os.path as op
import yaml
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression
import time
import pandas as pd
from joblib import Parallel, delayed


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
    reg_gene, target_gene = edge_gg.split('_')
    
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
        'source_gene': reg_gene, 'target_gene': target_gene, 'best_tx': best_tx,
        'S1_mean': S1_mean, 'S2_mean': S2_mean,
        'S1_median': S1_median, 'S2_median': S2_median,
        'E1': E1, 'E2': E2,
        'ratio': S2_mean / (S1_mean + eps) if E1 else np.inf,
        'source_category': category
    }


def canonical_vs_source_as(net1_tf,net2_tf):
    """
    The method iterates over all edges and computes if for a given gene-gene interaction,
    it is more likely that the source gene (as a whole) is the regulator, or if a specific isoform
    has better explainability for a given target.

    Args:
        net1_tf (pd.DataFrame): Gene-gene network
        net2_tf (pd.DataFrame): Tx-gene netork

    Returns:
        pd.DataFrame: Edge categorization of edge into source-isoform specific, source-gene specific, etc.
    """
    

    net1_tf['edge_gg'] = net1_tf['source_gene'] + '_' + net1_tf['target_gene']
    net2_tf['edge_gg'] = net2_tf['source_gene'] + '_' + net2_tf['target_gene']

    # Use MEAN importance for ratio calculation, MEDIAN for ranking
    # Net1: aggregate to gene-gene level
    # Can there even be duplicates: Imho no.
    net1_agg = net1_tf.groupby('edge_gg').agg({'mean_importance': 'max','median_importance': 'max'}).reset_index()
    net1_mean_imp = dict(zip(net1_agg['edge_gg'], net1_agg['mean_importance']))
    net1_median_imp = dict(zip(net1_agg['edge_gg'], net1_agg['median_importance']))

    # Net2: aggregate to gene-gene level (max across transcripts)
    # Compare against max edge: ok.
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
    sample_cols = [c for c in df.columns if c not in id_cols and pd.api.types.is_numeric_dtype(df[c])]
    
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
    
    # intialize empty
    result = {'tf_rho': np.nan, 'tf_pval': np.nan, 'tf_rho_conditional': np.nan, 'tf_pval_conditional': np.nan,
        'sf_rho': np.nan, 'sf_pval': np.nan, 'delta_usage': np.nan,'qc_ok': False, 'n_samples': 0}
    
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
        print('AAAA')
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




def tf_sf_disambigouation_fully_as_aware(net3_tfsf, regulator_list, transcript_data, usage_df, reliability_df,
                                         RHO_MIN = 0.3,  Q_MIN = 0.05, DU_MIN = 0.1, n_cores = 16):
    """
    Split the part of net3 into sets where there is evidence for the regulator acting as a splice factor,
    and where there is evidence for it being a transcription factor.

    Args:
        net3_tfsf (_type_): Net to be split into tf-like, sf-like etc.
        regulator_list (_type_): Regulator categories
        transcript_data (_type_): Transcript data (Needed to regress out TF data)
        usage_df (_type_): Fraction of counts attributed to this transcript
        reliability_df (_type_): Is transcript reliable.
        RHO_MIN (float, optional): _description_. Defaults to 0.3.
        Q_MIN (float, optional): _description_. Defaults to 0.05.
        DU_MIN (float, optional): _description_. Defaults to 0.1.

    Returns:
        _type_: _description_
    """
    # Step 1 get the data containing the sf and tf
    if len(net3_tfsf) > 0:
        sf_tx_list = regulator_list[regulator_list['Regulator_type'] == 'SF']['Transcript stable ID'].tolist()
        sf_tx_in_expr = [tx for tx in sf_tx_list if tx in transcript_data.index]
        if len(sf_tx_in_expr) > 0:
            sf_expr_matrix = transcript_data.loc[sf_tx_in_expr].T
            sf_expr_matrix = sf_expr_matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
        else:
            sf_expr_matrix = None
        
        tfsf_valid = net3_tfsf[net3_tfsf['source_transcript'].isin(transcript_data.index) & net3_tfsf['target_transcript'].isin(transcript_data.index)].copy()


    if len(tfsf_valid) > 0:
        tfsf_records = tfsf_valid.to_dict('records')
    
        # 2. Define the worker function for a chunk of rows
        def process_tfsf_chunk(chunk):
            results = []
            for row in chunk:
                # Compute the evidence metrics
                evidence = compute_tfsf_evidence(
                    row['source_transcript'], row['target_transcript'],
                    transcript_data, usage_df, reliability_df,
                    sf_expr_matrix
                )
                
                # Combine the original row data with the new evidence
                # This uses dictionary unpacking (**) for a clean merge
                results.append({
                    'source_transcript': row['source_transcript'],
                    'source_gene': row['source_gene'],
                    'target_transcript': row['target_transcript'],
                    'target_gene': row['target_gene'],
                    'mean_importance': row['mean_importance'],
                    'median_importance': row['median_importance'],
                    **evidence
                })
            return results

        chunks = np.array_split(tfsf_records, n_cores)

        # 4. Execute in parallel
        # n_jobs=16 distributes the chunks across your CPU cores
        parallel_results = Parallel(n_jobs=n_cores)(
            delayed(process_tfsf_chunk)(chunk) for chunk in chunks
        )

        # 5. Flatten the results and create the final DataFrame
        flattened_results = [item for sublist in parallel_results for item in sublist]
        set_d_full = pd.DataFrame(flattened_results)
        
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
        set_d_full.columns = ['source_transcript', 'source_gene', 'target_transcript', 'target_gene', 'mean_importance', 'median_importance']
        set_d_full['tfsf_category'] = 'tfsf_ambiguous'
    
    return set_d_full


        
def as_source_vs_as_full(net2_for_t2, net3_tf_only, r_iso=1.5, r_gene=1.5, r_eq=1.2, eps=1e-6):

    net2_for_t2['edge_tg'] = net2_for_t2['source_transcript'] + '_' + net2_for_t2['target_gene']
    net2_tg_mean_imp = dict(zip(net2_for_t2['edge_tg'], net2_for_t2['mean_importance']))

    net3_tf_only['edge_tg'] = net3_tf_only['source_transcript'] + '_' + net3_tf_only['target_gene']
    net3_tf_only['edge_tt'] = net3_tf_only['source_transcript'] + '_' + net3_tf_only['target_transcript']
    #Make mapping between the edges
    map_transcripts = dict(zip(net3_tf_only['edge_tt'], net3_tf_only['edge_tg']))

    net3_tg_mean = dict(zip(net3_tf_only['edge_tt'], net3_tf_only['mean_importance']))
    
    # Categorize TF edges (from Net2 TF+TF_SF vs Net3 TF)
    set_b_tf_rows = [categorize_target_resolution(et, map_transcripts[et], net2_tg_mean_imp, net3_tg_mean, r_iso=r_iso, r_gene=r_gene, r_eq=r_eq, eps=eps) for et in map_transcripts.keys()]
    set_b = pd.DataFrame(set_b_tf_rows)

    set_b['max_mean'] = set_b[['S2_mean', 'S3_mean']].max(axis=1)
    set_b = set_b.sort_values('max_mean', ascending=False)

    print(f"\nFinal Set B: {len(set_b):,} rows")
    #print(f"  - TF edges (from Net2 TF+TF_SF vs Net3 TF): {len(set_b_tf):,}")
    #print(f"  - TF_SF (tf_like from Set D): {len(set_b_tfsf) if len(set_b_tfsf) > 0 else 0:,}")
    print("\nCategory distribution:")
    for cat, count in set_b['target_category'].value_counts().items():
        print(f"  {cat:30} {count:>8,} ({count/len(set_b)*100:>5.1f}%)")
    print("\nBy regulator type:")

    return set_b   



def categorize_target_resolution(edge_tt, edge_tg, net2_mean, net3_mean, r_iso=1.5, r_gene=1.5, r_eq=1.2, eps=1e-6):
    """
    Categorize edge based on target resolution.
    Uses MEAN importance for ratio calculation (AlterNet 1.0 style).
    """
    # Get gene edge 
    S2_mean = net2_mean.get(edge_tg, 0)
    S3_mean = net3_mean.get(edge_tt, 0)
    
    E2, E3 = S2_mean > 0, S3_mean > 0
    reg_tx, target_gene = edge_tg.split('_')
    reg_tx, target_transcript = edge_tt.split('_')

    # Use MEAN importance for fold-change calculation
    if E3 and not E2:
        category = 'target_isoform_specific'
    elif E2 and not E3:
        category = 'target_gene_specific'
    elif E2 and E3:
        ratio = S3_mean / (S2_mean + eps)
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
        'source_transcript': reg_tx, 'target_transcript': target_transcript, 'target_gene': target_gene,
        'S2_mean': S2_mean, 'S3_mean': S3_mean,
        'E2': E2, 'E3': E3,
        'ratio': S3_mean / (S2_mean + eps) if E2 else np.inf,
        'target_category': category
    }




def get_presentation(row):
    """Derive target_tx_resolved, target_tx_dominant, n_target_tx from target_tx_set."""
    tx_set_str = row.get('target_tx_set', '')
    if pd.isna(tx_set_str) or tx_set_str == '':
        return pd.Series({'n_target_tx': 0, 'target_tx_resolved': '', 'target_tx_dominant': ''})
    
    tx_list = [t.strip() for t in str(tx_set_str).split(',') if t.strip()]
    n = len(tx_list)
    
    # Resolved: if exactly 1 transcript
    resolved = tx_list[0] if n == 1 else ''
    
    return pd.Series({'n_target_tx': n, 'target_tx_resolved': resolved, 'target_tx_dominant': ''})



def compute_sf_splicing_evidence(sf_tx, target_gene, target_tx_set, sf_vals, gene_to_all_transcripts, usage_df_indexed,
                                  transcript_mean_expr, reliability_df_indexed, epsilon=1e-6, q_min = 0.05, du_min=0.1, dom_min = 0.5, rho_min=0.3,
                                  min_isoforms_for_splicing = 2, top_m_expressed = 3):

    #Get all isoforms of the target gene
    all_gene_txs = gene_to_all_transcripts.get(target_gene, [])
    # Check of they are in the usage df and count available transcripts.
    all_gene_txs_in_data = [tx for tx in all_gene_txs if tx in usage_df_indexed.index]
    
    if len(all_gene_txs_in_data) < min_isoforms_for_splicing:
        return None
    # Check of the transcripts we passed are in usage df
    txs_to_test = set(tx for tx in target_tx_set if tx in usage_df_indexed.index)

    # Mean expression over transcripts in data for ranking.

    gene_txs_sorted = sorted(all_gene_txs_in_data, key=lambda tx: transcript_mean_expr.get(tx, 0), reverse=True)
    for tx in gene_txs_sorted[:top_m_expressed]:
        txs_to_test.add(tx)
    txs_to_test = list(txs_to_test)
    
    if len(txs_to_test) == 0:
        return None
    
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
        
        correlations.append({'source_transcript': sf_tx ,'target_transcript': tx, 'rho': rho, 'pval': pval, 'delta_usage': delta_usage})

    if len(correlations) == 0:
        return None
    
    _, qvals, _, _ = multipletests([c['pval'] for c in correlations], method='fdr_bh')

    correlations = pd.DataFrame(correlations)
    correlations['usage_reliable'] = True
    correlations['qval'] = qvals

    correlations['significant'] = (correlations['qval'] <= q_min) & (abs(correlations['rho']) >= rho_min) &  (abs(correlations['delta_usage']) >= du_min)
    n_sig = np.sum(correlations.significant)
    correlations['n_sig'] = n_sig
    
    if n_sig >= 2:
        # Push pull pattern if at least one SF isoform positively and one negatively correlated with target.
        # and significant and above threshold
        correlations['push_pull'] = (correlations['rho'].max() > 0) & (correlations['rho'].min() < 0)
    
    if n_sig > 0:
        abs_rhos = correlations['rho'].abs()
        correlations['dominance'] = correlations['rho'] / (abs_rhos.sum() + epsilon)
        correlations.loc[~correlations['significant'], 'dominance'] = np.nan

    correlations['sf_category'] = 'sf_expression_associated'
    if n_sig >= 2:
        correlations['sf_category'] = 'sf_splicing_supported_broad'
        evid = (correlations['push_pull']) | (correlations['dominance'] >= dom_min)
        correlations.loc[evid, 'sf_category'] =  'sf_splicing_supported_specific'

    return correlations


def compute_set_c(net3_sf, transcript_data, gene_to_all_transcripts, usage_df_indexed, reliability_df_indexed, sample_cols, epsilon=1e-6, n_cores=16):
    # Set C and Set D: PSI/Usage thresholds

    # Aggregate SF edges to (SF_tx, target_gene) level
    sf_edges = net3_sf.groupby(['source_transcript', 'source_gene', 'target_gene']).agg({
        'mean_importance': ['sum', 'max', 'count'],
        'median_importance': ['sum', 'max'],
        'target_transcript': lambda x: set(x)
    }).reset_index()
    sf_edges.columns = ['source_transcript', 'source_gene', 'target_gene','mean_importance_sum', 'mean_importance_max', 
                        'n_target_tx', 'median_importance_sum', 'median_importance_max', 'target_tx_set']

    sf_transcripts = sf_edges['source_transcript'].unique()
    transcript_temp = transcript_data.reset_index()
    transcript_mean_expr = transcript_temp[sample_cols].mean(axis=1).to_dict()

    sf_expr = transcript_temp[transcript_temp['transcript_id'].isin(sf_transcripts)].copy()
    sf_expr = sf_expr.set_index('transcript_id')[sample_cols]
    sf_expr_dict = {tx: np.array(sf_expr.loc[tx].values, dtype=np.float64) for tx in sf_expr.index}


    if len(sf_expr_dict) > 0:
        print(sf_edges.shape[0])
        edge_records = sf_edges.to_dict('records')

        flattened_results = []
        i = 0
        for row in edge_records:
            result = compute_sf_splicing_evidence(
                    row['source_transcript'], row['target_gene'], row['target_tx_set'],
                    sf_expr_dict[row['source_transcript']], gene_to_all_transcripts,
                    usage_df_indexed, transcript_mean_expr, reliability_df_indexed
                )
            flattened_results.append(result)
            i = i+1
            if i %1000 == 0 : print(i)
        set_c = pd.concat(flattened_results)

        # 1. Convert to a list of dicts for faster iteration if sticking to rows
        # edge_records = sf_edges.to_dict('records')

        # def process_chunk(chunk):
        #     results = []
        #     for row in chunk:
        #         result = compute_sf_splicing_evidence(
        #             row['source_transcript'], row['target_gene'], row['target_tx_set'],
        #             sf_expr_dict[row['source_transcript']], gene_to_all_transcripts,
        #             usage_df_indexed, transcript_mean_expr, reliability_df_indexed
        #         )
        #         results.append(result)
        #     print('Chunk done')
                
        #     return results # Return a list of DataFrames/dicts instead of concating inside

        # chunks = np.array_split(edge_records, n_cores)

        # # Use require='sharedmem' or prefer='threads' if compute_sf_splicing_evidence 
        # # utilizes underlying C-libraries (numpy, scipy, samtools) which release the GIL.
        # parallel_results = Parallel(n_jobs=n_cores, prefer="threads")(
        #     delayed(process_chunk)(chunk) for chunk in chunks
        # )

        # flattened_dfs = [df for sublist in parallel_results for df in sublist]
        # set_c = pd.concat(flattened_dfs, ignore_index=True)

    else:
        set_c = sf_edges.copy()
        set_c['sf_category'] = 'sf_expression_associated'

    set_c = set_c.merge(net3_sf, left_on = ['source_transcript', 'target_transcript'], right_on = ['source_transcript', 'target_transcript'])

    # Sort by median importance (AlterNet 1.0 style)
    set_c = set_c.sort_values('mean_importance', ascending=False)

    print(f"\nFinal Set C: {len(set_c):,} rows")

    print("\nCategory distribution:")
    for cat, count in set_c['sf_category'].value_counts().items():
        print(f"  {cat:35} {count:>8,} ({count/len(set_c)*100:>5.1f}%)")
    #print("\nBy regulator type:")
    #print(set_c.groupby(['reg_type', 'sf_category']).size().unstack(fill_value=0))

    return set_c

