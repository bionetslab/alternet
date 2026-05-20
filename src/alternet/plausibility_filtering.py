
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

def filter_set_a(set_a, gene_dominance, gene_n_isoforms,
                  dom_threshold=0.9, fc_eq_max=2.0, eps=1e-6):
    """
    Apply plausibility filters to Set A (source resolution).
    Updated for Part 2 v2 column names.
    """
    df = set_a.copy()
    
    # Add dominance and isoform count
    # change regulator to source!
    df['reg_dominance'] = df['source_gene'].map(gene_dominance).fillna(0.5)
    df['reg_n_isoforms'] = df['source_gene'].map(gene_n_isoforms).fillna(1)
    
    # Initialize filter columns
    df['is_plausible'] = True
    df['filter_reasons'] = ''
    
    # Identify importance columns (v2 uses S1_mean/S2_mean for ratios)
    s1_col = 'S1_mean'
    s2_col = 'S2_mean'
    
    # Compute ratio for reference (using mean importance for ratio - AlterNet 1.0 style)
    df['ratio_S2_S1'] = df[s2_col] / (df[s1_col] + eps)
    
    # Filter 1: Regulator equivalence (single isoform)
    # There is a single isoform it should be found in both nets (should be equivalent)
    # It is not, hence it is implausible.
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
    # Similar reasoning: it is dominant, but claimed to be isoform specific.
    # If there is a regulation specific to this isoform w
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



def filter_set_b(set_b, gene_dominance, gene_n_isoforms,
                  dom_threshold=0.9, fc_eq_max=2.0, eps=1e-6):
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
    s3_col = 'S3_mean'
    
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


def filter_set_c(set_c, gene_n_isoforms, min_tx_per_gene=2):
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



def filter_set_d(set_d, gene_dominance, gene_n_isoforms, 
                  dom_threshold=0.9, min_tx=2):
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