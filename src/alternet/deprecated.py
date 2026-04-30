
def filter_aggregated(data, top_percentile=0.3, threshold_frequency=5, importance_column='median_importance', frequency_column = 'frequency'):
    '''
    Filter out the top `threshold_importance` percent of data based on the importance column 
    and retain only edges that appear at least `threshold_frequency` times.

    Parameters:

        data : pd.DataFrame
            DataFrame containing the edges of the GRN inference, including importance values 
            and frequency of appearance after aggregation.

        threshold_importance : float
            Percentage (between 0 and 1) of top entries to retain based on the importance column.

        threshold_frequency : int
            Minimum number of times an edge must appear in the GRN inference to be retained.

        importance_column : str
            Name of the column containing importance values.

    Returns:
    
        pd.DataFrame
            Filtered DataFrame based on the specified importance and frequency thresholds.
    '''

    n_before = data.shape[0]
    #only select edges that have been found in threshold_frequency number of times in grn inference
    freq_mask = data[frequency_column] >= threshold_frequency
    
    # get value for threshold importance
    importance_threshold = data.loc[freq_mask, importance_column].quantile(1 - top_percentile)
    

    data =  data.loc[freq_mask & (data[importance_column] >= importance_threshold)]
    n_after = data.shape[0]

    filter_info = {'importance_freq': {'before': n_before, 'n_after': n_after}}
    return data.copy(), importance_threshold, filter_info