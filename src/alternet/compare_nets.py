import numpy as np
import pandas as pd



def filter_edges(df, max_p_value=0.05, importance_percentile=0.7,freq_col='pvalue_westfall_young', imp_col='importance'):
    """
    Filter edges by frequency and importance (AlterNet 1.0 style).
    Uses MEDIAN importance for filtering (more robust to outliers).
    """
    n_before = len(df)
    
    # Frequency filter first
    freq_mask = df[freq_col] < max_p_value
    df_filtered = df[freq_mask].copy()
    n_after_freq = len(df_filtered)
    
    # Importance filter (top percentile) - using median_importance
    importance_threshold = df_filtered[imp_col].quantile(importance_percentile)
    df_filtered = df_filtered[df_filtered[imp_col] >= importance_threshold].copy()
    n_after_imp = len(df_filtered)

    return df_filtered



def plausibility_filtering(df,gene_dominance, r_dom = 0.9):
    df['reg_dominance'] = df['source_gene'].map(gene_dominance).fillna(0.5)
    df['target_dominance'] = df['target_gene'].map(gene_dominance).fillna(0.5)

        # Initialize filter columns
    df['is_plausible'] = True

    # If both isoforms make a large proportion of the total expression, we should
    # at least find three edges for all more or less identical interactions.
    mask_iso_specific = ~((df.reg_dominance>r_dom) & (df.target_dominance>r_dom) & df['is_unique_edge'])
    df['is_plausible'] = mask_iso_specific

    return df



def canonical_names(df, tx2gene, gene2regtype):
    df = df.rename(columns = {'TF': 'source'})

    df["source_type"] = np.where(df["source"].str.startswith("ENSG"), "gene", "transcript")

    df["target_type"] = np.where(df["target"].str.startswith("ENSG"), "gene", "transcript")

    df["source_gene"] = np.where(df["source_type"] == 'gene', df["source"], df["source"].map(tx2gene))

    df["target_gene"] = np.where(df["target_type"] == 'gene', df["target"], df["target"].map(tx2gene))

    df["source_transcript"] = np.where(df["source_type"] == "transcript", df["source"], np.nan )
    
    df["target_transcript"] = np.where(df["target_type"] == "transcript", df["target"], np.nan)
    
    df["reg_type"] = df['source'].map(gene2regtype)


    return df


def get_best_variable(df, r_iso=0.66,r_eq=0.8, r_iso_strict=0.5, eps=1e-6):
    df['edge_gg'] = df['source_gene'] + '_' + df['target_gene']
    df['edge_type'] = df['source_type'] + '-' + df['target_type']

    # 1. Find the maximum 'importance' per group and map it to all elements
    df['max_importance'] = df.groupby('edge_gg')['importance'].transform('max')
    
    df["is_unique_edge"] = df.groupby("edge_gg")["edge_gg"].transform("count") == 1

    # 2. Compute the ratio of each element's importance to its group's maximum
    df['importance_ratio'] = df['importance'] / df['max_importance'] + eps

    is_max = df['importance'] == df['max_importance']

    df['max_other_ratio'] = (df['importance_ratio'].where(~is_max)
    .groupby(df['edge_gg'])
    .transform('max')
    .fillna(0)  # If a group only has 1 element, the max other ratio is 0
    )

        # 5. Define the conditions for each category
    # Specific: The max element is way ahead of the runner-up
    cond_specific_strict = is_max & (df['max_other_ratio'] < r_iso_strict)

    cond_specific = is_max & (df['max_other_ratio'] < r_iso)

    # Equivalent: There is no singular max because the runner-up is within the 1/r_eq window
    # This labels both the max and any other elements within that top window as 'equivalent'
    cond_equivalent = (df['importance_ratio'] >= r_eq) & (df['max_other_ratio'] >= r_eq)

    # 6. Apply categories using np.select
    conditions = [cond_specific_strict, cond_specific, cond_equivalent ]
    choices = ['specific_strict', 'specific', 'equivalent']

    df['category'] = np.select(conditions, choices, default='other')

    df['n_equivalent_edge_types'] = (df['edge_type']
        .where(cond_equivalent)
        .groupby(df['edge_gg'])
        .transform('nunique')
    )
    # Clean up helper columns
    df = df.drop(columns=['max_importance', 'max_other_ratio'])
    return df
