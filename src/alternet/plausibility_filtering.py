"""Plausibility filtering and annotation pipeline for AlterNet 2.0."""

import numpy as np
import pandas as pd


def compute_dominance_metrics(
    transcript_tpm: pd.DataFrame,
    sample_cols: list,
    min_tpm: float = 0.1,
    eps: float = 1e-6,
) -> tuple:
    """Compute dominance metrics for each gene.

    Args:
        transcript_tpm: Transcript TPM DataFrame with ``transcript_id``
            and ``gene_id`` columns.
        sample_cols: List of sample column names.
        min_tpm: Minimum mean TPM to consider a transcript expressed.
        eps: Small constant to avoid division by zero.

    Returns:
        Tuple of (gene_dominance, gene_n_isoforms, tx_expression_share) dicts.
    """
    tx_mean_tpm = transcript_tpm[sample_cols].mean(axis=1)
    df = transcript_tpm.copy()
    df["mean_tpm"] = tx_mean_tpm

    expressed = df[df["mean_tpm"] >= min_tpm].copy()
    gene_totals = expressed.groupby("gene_id")["mean_tpm"].sum()
    expressed["gene_total"] = expressed["gene_id"].map(gene_totals)
    expressed["expr_share"] = expressed["mean_tpm"] / (expressed["gene_total"] + eps)

    gene_dominance = expressed.groupby("gene_id")["expr_share"].max().to_dict()
    gene_n_isoforms = expressed.groupby("gene_id").size().to_dict()
    tx_expression_share = dict(zip(expressed["transcript_id"], expressed["expr_share"]))

    return gene_dominance, gene_n_isoforms, tx_expression_share


def filter_set_a(
    set_a: pd.DataFrame,
    gene_dominance: dict,
    gene_n_isoforms: dict,
    dom_threshold: float = 0.9,
    fc_eq_max: float = 2.0,
    eps: float = 1e-6,
) -> pd.DataFrame:
    """Apply plausibility filters to Set A (source resolution).

    Args:
        set_a: Set A DataFrame from edge categorization.
        gene_dominance: Dict mapping gene_id to max isoform expression share.
        gene_n_isoforms: Dict mapping gene_id to number of expressed isoforms.
        dom_threshold: Dominance threshold above which a gene is considered dominant.
        fc_eq_max: Maximum allowed fold-change for equivalent classification.
        eps: Small constant to avoid division by zero.

    Returns:
        Filtered DataFrame with ``is_plausible`` and ``filter_reasons`` columns.
    """
    df = set_a.copy()
    df["reg_dominance"] = df["regulator_gene"].map(gene_dominance).fillna(0.5)
    df["reg_n_isoforms"] = df["regulator_gene"].map(gene_n_isoforms).fillna(1)
    df["is_plausible"] = True
    df["filter_reasons"] = ""

    s1_col, s2_col = "S1_mean", "S2_mean"
    df["ratio_S2_S1"] = df[s2_col] / (df[s1_col] + eps)

    mask_single_iso = df["reg_n_isoforms"] == 1
    mask_iso_specific = df["source_category"].isin(
        ["source_isoform_specific", "source_gene_specific"]
    )
    implausible_single = mask_single_iso & mask_iso_specific
    df.loc[implausible_single, "is_plausible"] = False
    df.loc[implausible_single, "filter_reasons"] += "single_isoform_reg;"

    mask_dominant = df["reg_dominance"] >= dom_threshold
    mask_isoform_specific = df["source_category"] == "source_isoform_specific"
    weak_evidence = mask_dominant & mask_isoform_specific & (df["ratio_S2_S1"] < 2.0)
    df.loc[weak_evidence, "is_plausible"] = False
    df.loc[weak_evidence, "filter_reasons"] += "dominant_reg_weak_evidence;"

    mask_equivalent = df["source_category"] == "source_equivalent"
    ratio_ok = (df["ratio_S2_S1"] >= 1 / fc_eq_max) & (df["ratio_S2_S1"] <= fc_eq_max)
    ratio_ok = ratio_ok | (df[s1_col] == 0) | (df[s2_col] == 0)
    inconsistent_equiv = mask_equivalent & ~ratio_ok
    df.loc[inconsistent_equiv, "source_category"] = "source_ambiguous"
    df.loc[inconsistent_equiv, "filter_reasons"] += "equiv_ratio_out_of_bounds;"

    return df


def filter_set_b(
    set_b: pd.DataFrame,
    gene_dominance: dict,
    gene_n_isoforms: dict,
    dom_threshold: float = 0.9,
    fc_eq_max: float = 2.0,
    eps: float = 1e-6,
) -> pd.DataFrame:
    """Apply plausibility filters to Set B (target resolution).

    Args:
        set_b: Set B DataFrame from edge categorization.
        gene_dominance: Dict mapping gene_id to max isoform expression share.
        gene_n_isoforms: Dict mapping gene_id to number of expressed isoforms.
        dom_threshold: Dominance threshold above which a gene is considered dominant.
        fc_eq_max: Maximum allowed fold-change for equivalent classification.
        eps: Small constant to avoid division by zero.

    Returns:
        Filtered DataFrame with ``is_plausible`` and ``filter_reasons`` columns.
    """
    df = set_b.copy()
    df["tgt_dominance"] = df["target_gene"].map(gene_dominance).fillna(0.5)
    df["tgt_n_isoforms"] = df["target_gene"].map(gene_n_isoforms).fillna(1)
    df["is_plausible"] = True
    df["filter_reasons"] = ""

    s2_col, s3_col = "S2_mean", "S3_mean_sum"
    df["ratio_S3_S2"] = df[s3_col] / (df[s2_col] + eps)

    mask_single_iso = df["tgt_n_isoforms"] == 1
    mask_isoform_specific = df["target_category"] == "target_isoform_specific"
    implausible_single = mask_single_iso & mask_isoform_specific
    df.loc[implausible_single, "is_plausible"] = False
    df.loc[implausible_single, "filter_reasons"] += "single_isoform_tgt;"

    mask_dominant = df["tgt_dominance"] >= dom_threshold
    dominant_isoform_specific = mask_dominant & mask_isoform_specific
    df.loc[dominant_isoform_specific, "filter_reasons"] += "dominant_tgt_review;"

    mask_equivalent = df["target_category"] == "target_equivalent"
    ratio_ok = (df["ratio_S3_S2"] >= 1 / fc_eq_max) & (df["ratio_S3_S2"] <= fc_eq_max)
    ratio_ok = ratio_ok | (df[s2_col] == 0) | (df[s3_col] == 0)
    inconsistent_equiv = mask_equivalent & ~ratio_ok
    df.loc[inconsistent_equiv, "target_category"] = "target_ambiguous"
    df.loc[inconsistent_equiv, "filter_reasons"] += "equiv_ratio_out_of_bounds;"

    if "reg_type" in df.columns:
        tfsf_edges = df["reg_type"] == "TF_SF"
        df.loc[tfsf_edges, "filter_reasons"] += "tfsf_tf_like;"

    return df


def filter_set_c(
    set_c: pd.DataFrame,
    gene_n_isoforms: dict,
    min_tx_per_gene: int = 2,
) -> pd.DataFrame:
    """Apply plausibility filters to Set C (SF splicing evidence).

    Args:
        set_c: Set C DataFrame from edge categorization.
        gene_n_isoforms: Dict mapping gene_id to number of expressed isoforms.
        min_tx_per_gene: Minimum isoforms required for splicing claims.

    Returns:
        Filtered DataFrame with ``is_plausible`` and ``filter_reasons`` columns.
    """
    df = set_c.copy()
    df["tgt_n_isoforms"] = df["target_gene"].map(gene_n_isoforms).fillna(1)
    df["is_plausible"] = True
    df["filter_reasons"] = ""

    mask_splicing = df["sf_category"].str.contains("splicing_supported", na=False)
    mask_single_iso = df["tgt_n_isoforms"] < min_tx_per_gene
    invalid_splicing = mask_splicing & mask_single_iso
    df.loc[invalid_splicing, "sf_category"] = "sf_expression_associated"
    df.loc[invalid_splicing, "filter_reasons"] += "insufficient_isoforms_for_splicing;"

    mask_specific = df["sf_category"] == "sf_splicing_supported_specific"
    has_push_pull = df.get("push_pull", pd.Series(False, index=df.index)) == True
    has_multi_sig = df.get("n_sig", pd.Series(0, index=df.index)) >= 2
    invalid_specific = mask_specific & ~has_push_pull & ~has_multi_sig
    df.loc[invalid_specific, "sf_category"] = "sf_expression_associated"
    df.loc[invalid_specific, "filter_reasons"] += "specific_criteria_not_met;"

    unreliable = df.get("usage_reliable", pd.Series(True, index=df.index)) == False
    splicing_unreliable = mask_splicing & unreliable
    df.loc[splicing_unreliable, "sf_category"] = "sf_ambiguous"
    df.loc[splicing_unreliable, "filter_reasons"] += "usage_unreliable;"

    if "reg_type" in df.columns:
        tfsf_edges = df["reg_type"] == "TF_SF"
        df.loc[tfsf_edges, "filter_reasons"] += "tfsf_sf_like;"

    return df


def filter_set_d(
    set_d: pd.DataFrame,
    gene_dominance: dict,
    gene_n_isoforms: dict,
    dom_threshold: float = 0.9,
    min_tx: int = 2,
) -> pd.DataFrame:
    """Apply plausibility filters to Set D (TF+SF joint/ambiguous).

    Args:
        set_d: Set D DataFrame from edge categorization.
        gene_dominance: Dict mapping gene_id to max isoform expression share.
        gene_n_isoforms: Dict mapping gene_id to number of expressed isoforms.
        dom_threshold: Dominance threshold for target gene.
        min_tx: Minimum isoforms required for joint claims.

    Returns:
        Filtered DataFrame with ``is_plausible`` and ``filter_reasons`` columns.
    """
    df = set_d.copy()
    df["tgt_dominance"] = df["target_gene"].map(gene_dominance).fillna(0.5)
    df["tgt_n_isoforms"] = df["target_gene"].map(gene_n_isoforms).fillna(1)
    df["is_plausible"] = True
    df["filter_reasons"] = ""

    qc_failed = df.get("qc_ok", pd.Series(True, index=df.index)) == False
    df.loc[qc_failed, "tfsf_category"] = "tfsf_ambiguous"
    df.loc[qc_failed, "filter_reasons"] += "qc_failed;"

    mask_joint = df["tfsf_category"] == "tfsf_joint"
    mask_single_iso = df["tgt_n_isoforms"] < min_tx
    invalid_joint = mask_joint & mask_single_iso
    df.loc[invalid_joint, "tfsf_category"] = "tfsf_ambiguous"
    df.loc[invalid_joint, "filter_reasons"] += "joint_single_isoform;"

    mask_dominant_tgt = df["tgt_dominance"] >= dom_threshold
    dominant_joint = mask_dominant_tgt & mask_joint
    df.loc[dominant_joint, "filter_reasons"] += "dominant_tgt_joint_claim;"

    return df


def get_principal_tx(
    gene_id: str,
    digger_gene2tx: dict,
    tx_appris_label: dict,
) -> str:
    """Return the principal isoform for a gene using APPRIS labels.

    Args:
        gene_id: Ensembl gene ID.
        digger_gene2tx: Dict mapping gene_id to set of transcript IDs.
        tx_appris_label: Dict mapping transcript_id to APPRIS annotation label.

    Returns:
        Principal transcript ID, or ``None`` if not found.
    """
    candidates = digger_gene2tx.get(gene_id, set())
    best_tx, best_rank = None, 999
    for tx in candidates:
        label = tx_appris_label.get(tx, "")
        if "PRINCIPAL" in str(label).upper():
            parts = str(label).split(":")
            rank = int(parts[-1]) if len(parts) > 1 and parts[-1].isdigit() else 99
            if rank < best_rank:
                best_rank = rank
                best_tx = tx
    return best_tx


def get_unique_features(
    tx_id: str,
    gene_id: str,
    tx_exon_set: dict,
    tx_domain_set: dict,
    digger_gene2tx: dict,
    tx_appris_label: dict,
) -> tuple:
    """Return unique exons and domains for a transcript compared to baseline.

    Args:
        tx_id: Transcript ID.
        gene_id: Gene ID.
        tx_exon_set: Dict mapping transcript_id to set of exon ranks.
        tx_domain_set: Dict mapping transcript_id to set of Pfam domain IDs.
        digger_gene2tx: Dict mapping gene_id to set of transcript IDs.
        tx_appris_label: Dict mapping transcript_id to APPRIS annotation label.

    Returns:
        Tuple of (unique_exons_str, unique_domains_str).
    """
    if pd.isna(tx_id) or tx_id == "":
        return "", ""

    tx_exons = tx_exon_set.get(tx_id, set())
    tx_doms = tx_domain_set.get(tx_id, set())
    if not tx_exons and not tx_doms:
        return "", ""

    gene_txs = digger_gene2tx.get(gene_id, set()) if gene_id else set()
    principal = get_principal_tx(gene_id, digger_gene2tx, tx_appris_label)

    if principal and principal != tx_id:
        baseline_exons = tx_exon_set.get(principal, set())
        baseline_doms = tx_domain_set.get(principal, set())
    else:
        others = gene_txs - {tx_id}
        baseline_exons = (
            set().union(*(tx_exon_set.get(t, set()) for t in others))
            if others
            else set()
        )
        baseline_doms = (
            set().union(*(tx_domain_set.get(t, set()) for t in others))
            if others
            else set()
        )

    unique_exons = tx_exons - baseline_exons
    unique_doms = tx_doms - baseline_doms

    ue_str = ",".join(str(e) for e in sorted(unique_exons)) if unique_exons else ""
    ud_str = ",".join(sorted(unique_doms)) if unique_doms else ""
    return ue_str, ud_str


def get_shared_features(
    gene_id: str,
    digger_gene2tx: dict,
    tx_exon_set: dict,
    tx_domain_set: dict,
) -> tuple:
    """Return exons and domains shared by all isoforms of a gene.

    Args:
        gene_id: Ensembl gene ID.
        digger_gene2tx: Dict mapping gene_id to set of transcript IDs.
        tx_exon_set: Dict mapping transcript_id to set of exon ranks.
        tx_domain_set: Dict mapping transcript_id to set of Pfam domain IDs.

    Returns:
        Tuple of (shared_exons_str, shared_domains_str).
    """
    if pd.isna(gene_id) or gene_id == "":
        return "", ""

    gene_txs = digger_gene2tx.get(gene_id, set())
    if not gene_txs:
        return "", ""

    txs_with_data = [t for t in gene_txs if t in tx_exon_set]
    if not txs_with_data:
        return "", ""

    shared_exons = tx_exon_set.get(txs_with_data[0], set()).copy()
    shared_doms = tx_domain_set.get(txs_with_data[0], set()).copy()
    for t in txs_with_data[1:]:
        shared_exons &= tx_exon_set.get(t, set())
        shared_doms &= tx_domain_set.get(t, set())

    se_str = (
        ",".join(str(e) for e in sorted(shared_exons)) if shared_exons else ""
    )
    sd_str = ",".join(sorted(shared_doms)) if shared_doms else ""
    return se_str, sd_str


def annotate_digger_raw(
    tx_id: str,
    tx_exon_count: dict,
    tx_domain_count: dict,
    tx_domains: dict,
) -> tuple:
    """Return basic DIGGER annotation for a single transcript.

    Args:
        tx_id: Transcript ID.
        tx_exon_count: Dict mapping transcript_id to exon count.
        tx_domain_count: Dict mapping transcript_id to domain count.
        tx_domains: Dict mapping transcript_id to comma-separated domain string.

    Returns:
        Tuple of (exon_count, domain_count, domains_str).
    """
    if pd.isna(tx_id) or tx_id == "":
        return "", "", ""
    return (
        tx_exon_count.get(tx_id, ""),
        tx_domain_count.get(tx_id, ""),
        tx_domains.get(tx_id, ""),
    )


def _run_pipeline(
    set_a: pd.DataFrame,
    set_b: pd.DataFrame,
    set_c: pd.DataFrame,
    set_d: pd.DataFrame,
    set_b_unpacked: pd.DataFrame,
    set_c_unpacked: pd.DataFrame,
    transcript_tpm: pd.DataFrame,
    sample_cols: list,
    biomart: pd.DataFrame,
    appris_df: pd.DataFrame,
    digger_df: pd.DataFrame,
    dom_reg: float = 0.9,
    dom_tgt: float = 0.9,
    fc_eq_max: float = 2.0,
    min_tx_per_gene: int = 2,
    eps: float = 1e-6,
) -> dict:
    """Run the full plausibility filtering and annotation pipeline.

    Args:
        set_a: Set A DataFrame from edge categorization.
        set_b: Set B DataFrame from edge categorization.
        set_c: Set C DataFrame from edge categorization.
        set_d: Set D DataFrame from edge categorization.
        set_b_unpacked: Set B unpacked transcript-level DataFrame.
        set_c_unpacked: Set C unpacked transcript-level DataFrame.
        transcript_tpm: Transcript TPM DataFrame.
        sample_cols: List of sample column names.
        biomart: BioMart reference DataFrame.
        appris_df: APPRIS annotation DataFrame.
        digger_df: DIGGER annotation DataFrame.
        dom_reg: Dominance threshold for regulator genes.
        dom_tgt: Dominance threshold for target genes.
        fc_eq_max: Maximum fold-change for equivalent classification.
        min_tx_per_gene: Minimum isoforms for splicing claims.
        eps: Small constant to avoid division by zero.

    Returns:
        Dictionary with keys: ``set_a_filtered``, ``set_b_filtered``,
        ``set_c_filtered``, ``set_d_filtered``, ``set_b_unpacked``,
        ``set_c_unpacked``.
    """
    gene_dominance, gene_n_isoforms, tx_expr_share = compute_dominance_metrics(
        transcript_tpm, sample_cols, eps=eps
    )

    set_a_filtered = filter_set_a(set_a, gene_dominance, gene_n_isoforms, dom_reg, fc_eq_max, eps)
    set_b_filtered = filter_set_b(set_b, gene_dominance, gene_n_isoforms, dom_tgt, fc_eq_max, eps)
    set_c_filtered = filter_set_c(set_c, gene_n_isoforms, min_tx_per_gene)
    set_d_filtered = filter_set_d(set_d, gene_dominance, gene_n_isoforms, dom_tgt, min_tx_per_gene)

    # Propagate plausibility to unpacked tables
    set_b_unpacked = set_b_unpacked.copy()
    t2_plaus_key = set_b_filtered[
        ["regulator_tx", "target_gene", "is_plausible", "target_category"]
    ].copy()
    t2_plaus_key["edge_tg"] = (
        t2_plaus_key["regulator_tx"] + "|" + t2_plaus_key["target_gene"]
    )
    t2_plaus_dict = dict(zip(t2_plaus_key["edge_tg"], t2_plaus_key["is_plausible"]))
    t2_cat_dict = dict(zip(t2_plaus_key["edge_tg"], t2_plaus_key["target_category"]))
    if "regulator_tx" in set_b_unpacked.columns:
        set_b_unpacked["edge_tg"] = (
            set_b_unpacked["regulator_tx"] + "|" + set_b_unpacked["target_gene"]
        )
        set_b_unpacked["is_plausible"] = (
            set_b_unpacked["edge_tg"].map(t2_plaus_dict).fillna(False)
        )
        set_b_unpacked["target_category"] = set_b_unpacked["edge_tg"].map(t2_cat_dict).fillna(
            set_b_unpacked.get("target_category", "")
        )
        set_b_unpacked = set_b_unpacked.drop(columns=["edge_tg"])

    set_c_unpacked = set_c_unpacked.copy()
    t3_plaus_key = set_c_filtered[
        ["sf_tx", "target_gene", "is_plausible", "sf_category"]
    ].copy()
    t3_plaus_key["edge_sg"] = t3_plaus_key["sf_tx"] + "|" + t3_plaus_key["target_gene"]
    t3_plaus_dict = dict(zip(t3_plaus_key["edge_sg"], t3_plaus_key["is_plausible"]))
    t3_cat_dict = dict(zip(t3_plaus_key["edge_sg"], t3_plaus_key["sf_category"]))
    if "sf_tx" in set_c_unpacked.columns:
        set_c_unpacked["edge_sg"] = (
            set_c_unpacked["sf_tx"] + "|" + set_c_unpacked["target_gene"]
        )
        set_c_unpacked["is_plausible"] = (
            set_c_unpacked["edge_sg"].map(t3_plaus_dict).fillna(False)
        )
        set_c_unpacked["sf_category"] = set_c_unpacked["edge_sg"].map(t3_cat_dict).fillna(
            set_c_unpacked.get("sf_category", "")
        )
        set_c_unpacked = set_c_unpacked.drop(columns=["edge_sg"])

    # Build DIGGER lookups
    digger_df = digger_df.copy()
    digger_df["transcript_id"] = digger_df["Transcript stable ID"].str.split(".").str[0]

    tx_exon_count = (
        digger_df.groupby("transcript_id")["Exon rank in transcript"].nunique().to_dict()
    )
    domain_counts = (
        digger_df[digger_df["Pfam ID"].notna()]
        .groupby("transcript_id")["Pfam ID"]
        .nunique()
    )
    tx_domain_count = domain_counts.to_dict()

    def _get_domains(group):
        domains = group["Pfam ID"].dropna().unique()
        return ",".join(sorted(domains)) if len(domains) > 0 else ""

    tx_domains = digger_df.groupby("transcript_id").apply(_get_domains).to_dict()

    tx_exon_set = (
        digger_df.groupby("transcript_id")["Exon rank in transcript"]
        .apply(lambda x: set(x.dropna().unique()))
        .to_dict()
    )
    tx_domain_set = (
        digger_df[digger_df["Pfam ID"].notna()]
        .groupby("transcript_id")["Pfam ID"]
        .apply(lambda x: set(x.unique()))
        .to_dict()
    )
    digger_gene2tx = (
        digger_df.groupby(digger_df["Gene stable ID"].str.split(".").str[0])["transcript_id"]
        .apply(set)
        .to_dict()
    )

    # Build APPRIS lookup
    appris_df = appris_df.copy()
    appris_df["transcript_id"] = appris_df["Transcript ID"].str.split(".").str[0]
    tx_appris_label = dict(zip(appris_df["transcript_id"], appris_df["APPRIS Annotation"]))

    # BioMart lookups
    tx2gene = dict(zip(biomart["Transcript stable ID"], biomart["Gene stable ID"]))
    gene2symbol = dict(zip(biomart["Gene stable ID"], biomart["Gene name"]))
    tx_biotype = {
        tx: biomart[biomart["Transcript stable ID"] == tx]["Gene type"].values[0]
        if len(biomart[biomart["Transcript stable ID"] == tx]) > 0
        else ""
        for tx in list(set(biomart["Transcript stable ID"]))[:100]
    }
    # Faster version
    gene_biotype = dict(zip(biomart["Gene stable ID"], biomart["Gene type"]))
    tx_biotype = {tx: gene_biotype.get(gene, "") for tx, gene in tx2gene.items()}
    tx_gene_symbol = {tx: gene2symbol.get(gene, "") for tx, gene in tx2gene.items()}

    def _annotate_digger_raw_fn(tx_id):
        return annotate_digger_raw(tx_id, tx_exon_count, tx_domain_count, tx_domains)

    def _get_unique_features_fn(tx_id, gene_id):
        return get_unique_features(
            tx_id, gene_id, tx_exon_set, tx_domain_set, digger_gene2tx, tx_appris_label
        )

    def _get_shared_features_fn(gene_id):
        return get_shared_features(gene_id, digger_gene2tx, tx_exon_set, tx_domain_set)

    # Annotate Set A
    set_a_filtered["reg_appris"] = set_a_filtered["best_tx"].map(
        lambda x: tx_appris_label.get(x, "") if pd.notna(x) else ""
    )
    set_a_filtered["reg_biotype"] = set_a_filtered["best_tx"].map(
        lambda x: tx_biotype.get(x, "") if pd.notna(x) else ""
    )
    set_a_filtered["reg_expr_share"] = set_a_filtered["best_tx"].map(
        lambda x: tx_expr_share.get(x, np.nan) if pd.notna(x) else np.nan
    )
    set_a_filtered["reg_symbol"] = set_a_filtered["regulator_gene"].map(gene2symbol)
    set_a_filtered["tgt_symbol"] = set_a_filtered["target_gene"].map(gene2symbol)

    raw = set_a_filtered["best_tx"].apply(_annotate_digger_raw_fn)
    set_a_filtered["reg_exon_count"] = raw.apply(lambda x: x[0])
    set_a_filtered["reg_domain_count"] = raw.apply(lambda x: x[1])
    set_a_filtered["reg_domains"] = raw.apply(lambda x: x[2])

    def _digger_set_a(row):
        cat = row.get("source_category", "")
        tx = row.get("best_tx", "")
        gene = row.get("regulator_gene", "")
        if cat == "source_isoform_specific":
            ue, ud = _get_unique_features_fn(tx, gene)
            return pd.Series(
                {"reg_unique_exons": ue, "reg_unique_domains": ud,
                 "reg_shared_exons": "", "reg_shared_domains": ""}
            )
        elif cat == "source_gene_specific":
            se, sd = _get_shared_features_fn(gene)
            return pd.Series(
                {"reg_unique_exons": "", "reg_unique_domains": "",
                 "reg_shared_exons": se, "reg_shared_domains": sd}
            )
        return pd.Series(
            {"reg_unique_exons": "", "reg_unique_domains": "",
             "reg_shared_exons": "", "reg_shared_domains": ""}
        )

    digger_cols = set_a_filtered.apply(_digger_set_a, axis=1)
    set_a_filtered = pd.concat([set_a_filtered, digger_cols], axis=1)

    # Annotate Set B
    set_b_filtered["reg_appris"] = set_b_filtered["regulator_tx"].map(
        lambda x: tx_appris_label.get(x, "") if pd.notna(x) else ""
    )
    set_b_filtered["reg_biotype"] = set_b_filtered["regulator_tx"].map(
        lambda x: tx_biotype.get(x, "") if pd.notna(x) else ""
    )
    set_b_filtered["reg_symbol"] = set_b_filtered["regulator_gene"].map(gene2symbol)
    set_b_filtered["tgt_symbol"] = set_b_filtered["target_gene"].map(gene2symbol)

    raw = set_b_filtered["regulator_tx"].apply(_annotate_digger_raw_fn)
    set_b_filtered["reg_exon_count"] = raw.apply(lambda x: x[0])
    set_b_filtered["reg_domain_count"] = raw.apply(lambda x: x[1])
    set_b_filtered["reg_domains"] = raw.apply(lambda x: x[2])

    def _digger_set_b(row):
        cat = row.get("target_category", "")
        tx = row.get("target_tx_resolved", "")
        gene = row.get("target_gene", "")
        if cat == "target_isoform_specific":
            ue, ud = _get_unique_features_fn(tx, gene)
            return pd.Series(
                {"tgt_unique_exons": ue, "tgt_unique_domains": ud,
                 "tgt_shared_exons": "", "tgt_shared_domains": ""}
            )
        elif cat == "target_gene_specific":
            se, sd = _get_shared_features_fn(gene)
            return pd.Series(
                {"tgt_unique_exons": "", "tgt_unique_domains": "",
                 "tgt_shared_exons": se, "tgt_shared_domains": sd}
            )
        return pd.Series(
            {"tgt_unique_exons": "", "tgt_unique_domains": "",
             "tgt_shared_exons": "", "tgt_shared_domains": ""}
        )

    digger_cols = set_b_filtered.apply(_digger_set_b, axis=1)
    set_b_filtered = pd.concat([set_b_filtered, digger_cols], axis=1)

    # Annotate Set C
    set_c_filtered["sf_appris"] = set_c_filtered["sf_tx"].map(
        lambda x: tx_appris_label.get(x, "") if pd.notna(x) else ""
    )
    set_c_filtered["sf_biotype"] = set_c_filtered["sf_tx"].map(
        lambda x: tx_biotype.get(x, "") if pd.notna(x) else ""
    )
    set_c_filtered["sf_symbol"] = set_c_filtered["sf_gene"].map(gene2symbol)
    set_c_filtered["tgt_symbol"] = set_c_filtered["target_gene"].map(gene2symbol)

    raw = set_c_filtered["sf_tx"].apply(_annotate_digger_raw_fn)
    set_c_filtered["sf_exon_count"] = raw.apply(lambda x: x[0])
    set_c_filtered["sf_domain_count"] = raw.apply(lambda x: x[1])
    set_c_filtered["sf_domains"] = raw.apply(lambda x: x[2])

    def _digger_set_c(row):
        cat = row.get("sf_category", "")
        tx = (
            row.get("best_tx", "")
            if pd.notna(row.get("best_tx", ""))
            else row.get("target_tx_resolved", "")
        )
        gene = row.get("target_gene", "")
        if cat == "sf_splicing_supported_specific":
            ue, ud = _get_unique_features_fn(tx, gene)
            return pd.Series(
                {"tgt_unique_exons": ue, "tgt_unique_domains": ud,
                 "tgt_shared_exons": "", "tgt_shared_domains": ""}
            )
        elif cat == "sf_splicing_supported_broad":
            se, sd = _get_shared_features_fn(gene)
            return pd.Series(
                {"tgt_unique_exons": "", "tgt_unique_domains": "",
                 "tgt_shared_exons": se, "tgt_shared_domains": sd}
            )
        return pd.Series(
            {"tgt_unique_exons": "", "tgt_unique_domains": "",
             "tgt_shared_exons": "", "tgt_shared_domains": ""}
        )

    digger_cols = set_c_filtered.apply(_digger_set_c, axis=1)
    set_c_filtered = pd.concat([set_c_filtered, digger_cols], axis=1)

    return {
        "set_a_filtered": set_a_filtered,
        "set_b_filtered": set_b_filtered,
        "set_c_filtered": set_c_filtered,
        "set_d_filtered": set_d_filtered,
        "set_b_unpacked": set_b_unpacked,
        "set_c_unpacked": set_c_unpacked,
    }
