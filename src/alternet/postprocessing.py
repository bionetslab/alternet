import pandas as pd 
import os.path as op
import copy
import numpy as np


#++
def filter_edges(df, min_frequency=10, importance_percentile=0.7, freq_col='frequency', imp_col='median_importance'):
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

    return df_filtered

 
#++_++
def plausibility_filtering(isoform_unique, isoform_categories):

    n_before = isoform_unique.shape[0]
    isoform_unique = isoform_unique.merge(isoform_categories, left_on='source_transcript', right_on='Transcript stable ID')

    
    # remove all other dominant edges
    isoform_unique = isoform_unique[~isoform_unique.isoform_category.isin(['single'])]
    n_after_single = isoform_unique.shape[0]

    isoform_unique = isoform_unique[~isoform_unique.isoform_category.isin(['dominant'])]
    n_after_dominant = isoform_unique.shape[0]

    isoform_unique = isoform_unique.sort_values('median_importance', ascending=False)

    filter_info = {'n_before': n_before, 'n_after_equivalence': n_after_single,  'n_after_dominance': n_after_dominant}

    return isoform_unique, filter_info

 