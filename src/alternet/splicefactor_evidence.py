
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

def compute_sf_splicing_evidence(sf_tx, target_gene, target_tx_set, sf_vals, gene_to_all_transcripts, usage_df_indexed,
                                  transcript_mean_expr, reliability_df_indexed, epsilon=1e-6, q_min = 0.05, du_min=0.1, dom_min = 0.5, rho_min=0.3,
   
                                  min_isoforms_for_splicing = 2, top_m_expressed = 3):


    empty_result = {'source_transcript': sf_tx, 'target_transcript': list(target_tx_set)[0], 'rho': np.nan, 'pval': np.nan, 'delta_usage': np.nan,
       'usage_reliable': False, 'qval': np.nan, 'significant': False, 'n_sig': 0, 'dominance': 0,
       'sf_category': 'expression_associated', 'push_pull': False}
    
    empty_result = pd.DataFrame(empty_result, index = [0])

    #Get all isoforms of the target gene
    all_gene_txs = gene_to_all_transcripts.get(target_gene, [])
    # Check of they are in the usage df and count available transcripts.
    all_gene_txs_in_data = [tx for tx in all_gene_txs if tx in usage_df_indexed.index]
    
    if len(all_gene_txs_in_data) < min_isoforms_for_splicing:
        return empty_result
    # Check of the transcripts we passed are in usage df
    txs_to_test = set(tx for tx in target_tx_set if tx in usage_df_indexed.index)

    # Mean expression over transcripts in data for ranking.

    gene_txs_sorted = sorted(all_gene_txs_in_data, key=lambda tx: transcript_mean_expr.get(tx, 0), reverse=True)
    for tx in gene_txs_sorted[:top_m_expressed]:
        txs_to_test.add(tx)
    txs_to_test = list(txs_to_test)
    
    if len(txs_to_test) == 0:

        return empty_result
    
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
        print('No corrletation')
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

    counter = 0
    if len(sf_expr_dict) > 0:
        edge_records = sf_edges.to_dict('records')
        flattened_results = []
        i = 0
        for row in edge_records:
            result = compute_sf_splicing_evidence(
                    row['source_transcript'], row['target_gene'], row['target_tx_set'],
                    sf_expr_dict[row['source_transcript']], gene_to_all_transcripts,
                    usage_df_indexed, transcript_mean_expr, reliability_df_indexed
                )
            if result is None:
                counter = counter+1

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
    print(set_c.columns)
    print(f"Failing trnascripts {counter}")

    print(f'Final set c size{set_c.shape[0]}')
    set_c = set_c.merge(net3_sf, left_on = ['source_transcript', 'target_transcript'], right_on = ['source_transcript', 'target_transcript'])


    if 'delta_usage_x' in set_c.columns:
        # Happens if SF_TF has been disambiguated before
        set_c = set_c.drop(columns = {'delta_usage_x'})
        set_c = set_c.rename(columns = {'delta_usage_y': 'delta_usage'})

    # # Sort by median importance (AlterNet 1.0 style)
    set_c = set_c.sort_values('mean_importance', ascending=False)

    # print(f"\nFinal Set C: {len(set_c):,} rows")

    print("\nCategory distribution:")
    for cat, count in set_c['sf_category'].value_counts().items():
        print(f"  {cat:35} {count:>8,} ({count/len(set_c)*100:>5.1f}%)")
    #print("\nBy regulator type:")
    #print(set_c.groupby(['reg_type', 'sf_category']).size().unstack(fill_value=0))

    return set_c



def tf_sf_disambigouation_fully_as_aware(net3_tfsf, regulator_list, transcript_data, usage_df, reliability_df,
                                         RHO_MIN = 0.3,  Q_MIN = 0.05, DU_MIN = 0.1, n_cores = 16, transcript_col = 'transcript_id'):
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
        sf_tx_in_expr = [tx for tx in sf_tx_list if tx in transcript_data.index.values]
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
                results.append(row|evidence)
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



def compute_tfsf_evidence(reg_tx, target_tx, expr_df_indexed, usage_df_indexed, 
                          reliability_df_indexed, sf_expr_matrix):
    
    # intialize empty
    result = {'tf_rho': np.nan, 'tf_pval': np.nan, 'tf_rho_conditional': np.nan, 'tf_pval_conditional': np.nan,
        'sf_rho': np.nan, 'sf_pval': np.nan, 'delta_usage': np.nan,'qc_ok': False, 'n_samples': 0}
    
    if reg_tx not in expr_df_indexed.index or target_tx not in expr_df_indexed.index:
        print('E')
        return result
    if target_tx not in usage_df_indexed.index:
        print('A')
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


def calculate_transcript_usage(transcript_data, gene_col = 'gene_id', transcript_col = 'transcript_id',  min_gene_tpm=1.0, epsilon=1e-8):
    """
    Calculates transcript usage (fraction of gene expression) and reliability masks.
    
    Returns:
        usage_df: DataFrame of transcript usage (0.0 to 1.0)
        reliability_df: Binary mask (1.0 if gene expression >= min_gene_tpm)
        metadata: Dictionary containing gene-to-transcript mapping and mean expressions
    """
    genes = transcript_data[gene_col].values
    sample_cols = [c for c in transcript_data.columns if c not in [transcript_col, gene_col]]

    tpm_numeric = transcript_data[sample_cols].astype(float)
    gene_numeric_total  = transcript_data[[gene_col]+sample_cols].groupby(gene_col).transform('sum').astype(float)

    usage_values = tpm_numeric / (gene_numeric_total + epsilon)
    usage_df = pd.concat([transcript_data[[transcript_col, gene_col]], usage_values], axis=1)

    reliability_df = (gene_numeric_total >= min_gene_tpm).astype(float)
    reliability_df.index = transcript_data[transcript_col].values

    return usage_df, reliability_df


def compute_dominance_metrics(transcript_tpm, sample_cols, min_tpm=0.1, epsilon = 1e-6):
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
    expressed['expr_share'] = expressed['mean_tpm'] / (expressed['gene_total'] + epsilon)
    
    # Dominance = max expression share per gene
    gene_dominance = expressed.groupby('gene_id')['expr_share'].max().to_dict()
    
    # Number of expressed isoforms per gene
    gene_n_isoforms = expressed.groupby('gene_id').size().to_dict()
    
    # Transcript expression share
    tx_expression_share = dict(zip(expressed['transcript_id'], expressed['expr_share']))
    
    return gene_dominance, gene_n_isoforms, tx_expression_share