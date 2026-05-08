import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np

def create_hybrid_data(transcript_data, gene_data, tf_list, biomart_column='Transcript stable ID'):
    hybrid_data = pd.concat([transcript_data.loc[:, transcript_data.columns.isin(tf_list['Transcript stable ID'])], gene_data], axis =1)
    return hybrid_data
    




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