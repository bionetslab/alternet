import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np

def create_hybrid_data(transcript_data, gene_data, tf_list, biomart_column='Transcript stable ID'):
    hybrid_data = pd.concat([transcript_data.loc[:, transcript_data.columns.isin(tf_list['Transcript stable ID'])], gene_data], axis =1)
    return hybrid_data
    

def create_hybrid_data(transcript_data, gene_data, tf_list, biomart_column='Transcript stable ID'):
    hybrid_data = pd.concat([transcript_data.loc[:, transcript_data.columns.isin(tf_list[biomart_column])], gene_data], axis =1)
    return hybrid_data
    

def remove_low_expession_transcripts(data):
    data =data.set_index(['transcript_id', 'gene_id'])
    threshold = data.shape[1] * 0.1
    mask = (data == 0).sum(axis=1) > threshold
    data = data[~mask]
    data = data.reset_index()
    return data.copy()

def standardize_dataframe(df):
    """
    Standardizes a DataFrame column-wise using sklearn's StandardScaler.
    
    Parameters:
        df (pd.DataFrame): The input DataFrame.
        
    Returns:
        pd.DataFrame: The standardized DataFrame.
    """
    scaler = StandardScaler()

    scaled_array = scaler.fit_transform(df)
    
    standardized_df = pd.DataFrame(
        scaled_array, 
        columns=df.columns,
        index=df.index
    )
    return standardized_df




def variance_filtering(transcript_data, variance_percentile = 0.7, non_sample_cols = ['transcript_id', 'gene_id']):
    sample_cols = [c for c in transcript_data.columns if c not in non_sample_cols]
    expression_values = transcript_data[sample_cols].values
    log_expr = np.log1p(expression_values)
    variances = np.var(log_expr, axis=1)
    variance_threshold = np.quantile(variances, variance_percentile)

    n_before = len(transcript_data)
    transcript_data = transcript_data[variances > variance_threshold].copy()
    #print(f"Variance filter: {n_before} → {len(transcript_data)} transcripts")
    return transcript_data


def select_supported_isoforms(transcript_data, min_tx_tpm=1.0, min_gene_tpm=1.0, min_usage=0.1,
                              min_sample_frac=0.25, min_usage_sd=0.0, min_supported_isoforms=1,
                              gene_col='gene_id', transcript_col='transcript_id'):
    """
    Select isoforms with support within a single condition (no group labels, no DTU).

    A transcript is selected if, in at least `min_sample_frac` of the samples,
      - its gene is expressed (gene TPM >= min_gene_tpm),
      - the transcript is expressed (TPM >= min_tx_tpm),
      - the transcript makes up a relevant share of its gene (usage >= min_usage),
    and its usage varies across samples with expressed gene (SD >= min_usage_sd; 0 disables).
    Genes with fewer than `min_supported_isoforms` passing isoforms lose all their isoforms.

    Usage is computed against the gene total of ALL transcripts, so gene-level
    expression should be built from the unfiltered data.

    Returns:
        selected: transcript_data restricted to selected transcripts
        report: per-transcript metrics and pass flags (indexed by transcript_id)
    """
    sample_cols = [c for c in transcript_data.columns if c not in [transcript_col, gene_col]]
    tpm = transcript_data.set_index(transcript_col)[sample_cols].astype(float)
    genes = transcript_data[gene_col].values

    gene_tpm = tpm.groupby(genes).transform('sum')
    gene_expressed = gene_tpm >= min_gene_tpm
    usage = (tpm / gene_tpm.where(gene_expressed)).where(gene_expressed)

    min_samples = max(1, int(np.ceil(min_sample_frac * len(sample_cols))))

    report = pd.DataFrame({
        gene_col: genes,
        'n_samples_gene_expr': gene_expressed.sum(axis=1).values,
        'n_samples_tx_expr': (tpm >= min_tx_tpm).sum(axis=1).values,
        'n_samples_usage': (usage >= min_usage).sum(axis=1).values,
        'median_tpm': tpm.median(axis=1).values,
        'mean_usage': usage.mean(axis=1).values,
        'usage_sd': usage.std(axis=1).values,
    }, index=tpm.index)

    report['pass_gene_expr'] = report['n_samples_gene_expr'] >= min_samples
    report['pass_tx_expr'] = report['n_samples_tx_expr'] >= min_samples
    report['pass_usage'] = report['n_samples_usage'] >= min_samples
    report['pass_usage_sd'] = report['usage_sd'].fillna(0) >= min_usage_sd
    supported = report[['pass_gene_expr', 'pass_tx_expr', 'pass_usage', 'pass_usage_sd']].all(axis=1)

    report['n_supported_isoforms_gene'] = supported.groupby(report[gene_col]).transform('sum')
    report['n_isoforms_gene'] = report.groupby(gene_col)[gene_col].transform('size')
    report['selected'] = supported & (report['n_supported_isoforms_gene'] >= min_supported_isoforms)

    selected = transcript_data[transcript_data[transcript_col].isin(report.index[report['selected']])].copy()
    return selected, report


def remove_problematic_transcripts(transcript_data_scaled, gene_data_scaled, transcript_data_matrix):
    # Remove problematic transcripts (NaN or zero variance)
    nan_cols = transcript_data_scaled.columns[transcript_data_scaled.isna().any()].tolist()
    zero_var_cols = transcript_data_matrix.columns[transcript_data_matrix.std() == 0].tolist()
    bad_transcripts = set(nan_cols + zero_var_cols)
    print(bad_transcripts)
    if bad_transcripts:
        good_transcripts = [c for c in transcript_data_scaled.columns if c not in bad_transcripts]
        transcript_data_scaled = transcript_data_scaled[good_transcripts]
        transcript_data_matrix = transcript_data_matrix[good_transcripts]
        print(f"Removed {len(bad_transcripts)} problematic transcripts")

    # Fill remaining NaN
    transcript_data_scaled = transcript_data_scaled.fillna(0)
    gene_data_scaled = gene_data_scaled.fillna(0)


    return transcript_data_scaled, gene_data_scaled, transcript_data_matrix