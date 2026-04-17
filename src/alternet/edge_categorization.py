"""Edge categorization pipeline for AlterNet 2.0."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests


def filter_edges(
    df: pd.DataFrame,
    min_frequency: int = 10,
    importance_percentile: float = 0.7,
    freq_col: str = "frequency",
    imp_col: str = "median_importance",
) -> pd.DataFrame:
    """Filter edges by frequency and importance thresholds.

    Args:
        df: Edge DataFrame.
        min_frequency: Minimum edge frequency to keep.
        importance_percentile: Quantile threshold for importance filtering.
        freq_col: Column name for frequency.
        imp_col: Column name for importance (used for percentile filtering).

    Returns:
        Filtered DataFrame.
    """
    freq_mask = df[freq_col] >= min_frequency
    df_filtered = df[freq_mask].copy()

    importance_threshold = df_filtered[imp_col].quantile(importance_percentile)
    df_filtered = df_filtered[df_filtered[imp_col] >= importance_threshold].copy()

    return df_filtered


def categorize_source_resolution(
    edge_gg: str,
    net1_mean: dict,
    net2_mean: dict,
    net1_median: dict,
    net2_median: dict,
    net2_best_tx: dict,
    r_iso: float = 1.5,
    r_gene: float = 1.5,
    r_eq: float = 1.2,
    eps: float = 1e-6,
) -> dict:
    """Categorize an edge based on source (regulator) resolution.

    Compares Network 1 (gene-level) vs Network 2 (transcript-level) importance
    to determine if the source exhibits gene-specific, isoform-specific,
    equivalent, or ambiguous regulation.

    Args:
        edge_gg: Edge key in format ``'regulator_gene|target_gene'``.
        net1_mean: Dict of edge -> mean importance in Network 1.
        net2_mean: Dict of edge -> mean importance in Network 2.
        net1_median: Dict of edge -> median importance in Network 1.
        net2_median: Dict of edge -> median importance in Network 2.
        net2_best_tx: Dict of edge -> best transcript in Network 2.
        r_iso: Fold-change threshold for isoform-specific classification.
        r_gene: Fold-change threshold for gene-specific classification.
        r_eq: Equivalence band threshold.
        eps: Small constant to avoid division by zero.

    Returns:
        Dictionary with categorization results.
    """
    S1_mean = net1_mean.get(edge_gg, 0)
    S2_mean = net2_mean.get(edge_gg, 0)
    S1_median = net1_median.get(edge_gg, 0)
    S2_median = net2_median.get(edge_gg, 0)
    best_tx = net2_best_tx.get(edge_gg, None)

    E1, E2 = S1_mean > 0, S2_mean > 0
    reg_gene, target_gene = edge_gg.split("|")

    if E2 and not E1:
        category = "source_isoform_specific"
    elif E1 and not E2:
        category = "source_gene_specific"
    elif E1 and E2:
        ratio = S2_mean / (S1_mean + eps)
        if ratio >= r_iso:
            category = "source_isoform_specific"
        elif ratio <= 1 / r_gene:
            category = "source_gene_specific"
        elif 1 / r_eq <= ratio <= r_eq:
            category = "source_equivalent"
        else:
            category = "source_ambiguous"
    else:
        category = "source_ambiguous"

    return {
        "regulator_gene": reg_gene,
        "target_gene": target_gene,
        "best_tx": best_tx,
        "S1_mean": S1_mean,
        "S2_mean": S2_mean,
        "S1_median": S1_median,
        "S2_median": S2_median,
        "E1": E1,
        "E2": E2,
        "ratio": S2_mean / (S1_mean + eps) if E1 else np.inf,
        "source_category": category,
    }


def categorize_target_resolution(
    edge_tg: str,
    net2_mean: dict,
    net3_mean: dict,
    net3_max: dict,
    net3_dom: dict,
    net2_median: dict,
    net3_median: dict,
    tx2gene: dict,
    r_iso: float = 1.5,
    r_gene: float = 1.5,
    r_eq: float = 1.2,
    eps: float = 1e-6,
) -> dict:
    """Categorize an edge based on target resolution.

    Compares Network 2 (gene-level targets) vs Network 3 (transcript-level
    targets) to determine target specificity.

    Args:
        edge_tg: Edge key in format ``'regulator_tx|target_gene'``.
        net2_mean: Dict of edge -> mean importance in Network 2.
        net3_mean: Dict of edge -> summed mean importance in Network 3.
        net3_max: Dict of edge -> max mean importance in Network 3.
        net3_dom: Dict of edge -> dominance score in Network 3.
        net2_median: Dict of edge -> median importance in Network 2.
        net3_median: Dict of edge -> summed median importance in Network 3.
        tx2gene: Dict mapping transcript ID to gene ID.
        r_iso: Fold-change threshold for isoform-specific classification.
        r_gene: Fold-change threshold for gene-specific classification.
        r_eq: Equivalence band threshold.
        eps: Small constant to avoid division by zero.

    Returns:
        Dictionary with categorization results.
    """
    S2_mean = net2_mean.get(edge_tg, 0)
    S3_mean_sum = net3_mean.get(edge_tg, 0)
    S3_mean_max = net3_max.get(edge_tg, 0)
    S2_median = net2_median.get(edge_tg, 0)
    S3_median = net3_median.get(edge_tg, 0)
    dominance = net3_dom.get(edge_tg, 0)

    E2, E3 = S2_mean > 0, S3_mean_sum > 0
    reg_tx, target_gene = edge_tg.split("|")
    reg_gene = tx2gene.get(reg_tx, None)

    if E3 and not E2:
        category = "target_isoform_specific"
    elif E2 and not E3:
        category = "target_gene_specific"
    elif E2 and E3:
        ratio = S3_mean_sum / (S2_mean + eps)
        if ratio >= r_iso:
            category = "target_isoform_specific"
        elif ratio <= 1 / r_gene:
            category = "target_gene_specific"
        elif 1 / r_eq <= ratio <= r_eq:
            category = "target_equivalent"
        else:
            category = "target_ambiguous"
    else:
        category = "target_ambiguous"

    return {
        "regulator_tx": reg_tx,
        "regulator_gene": reg_gene,
        "target_gene": target_gene,
        "S2_mean": S2_mean,
        "S3_mean_sum": S3_mean_sum,
        "S3_mean_max": S3_mean_max,
        "S2_median": S2_median,
        "S3_median": S3_median,
        "dominance": dominance,
        "E2": E2,
        "E3": E3,
        "ratio": S3_mean_sum / (S2_mean + eps) if E2 else np.inf,
        "target_category": category,
    }


def compute_tfsf_evidence(
    reg_tx: str,
    target_tx: str,
    expr_df_indexed: pd.DataFrame,
    usage_df_indexed: pd.DataFrame,
    reliability_df_indexed: pd.DataFrame,
    sf_expr_matrix,
    sample_cols: list,
    rho_min: float = 0.3,
    q_min: float = 0.05,
    du_min: float = 0.1,
) -> dict:
    """Compute TF/SF evidence for a TF_SF regulator - target transcript pair.

    Args:
        reg_tx: Regulator transcript ID.
        target_tx: Target transcript ID.
        expr_df_indexed: Transcript expression DataFrame indexed by transcript ID.
        usage_df_indexed: Isoform usage DataFrame indexed by transcript ID.
        reliability_df_indexed: Usage reliability DataFrame indexed by transcript ID.
        sf_expr_matrix: SF expression matrix (samples x SF transcripts) or None.
        sample_cols: List of sample column names.
        rho_min: Minimum Spearman correlation threshold.
        q_min: Maximum FDR-adjusted p-value threshold.
        du_min: Minimum delta usage threshold.

    Returns:
        Dictionary with evidence metrics.
    """
    result = {
        "tf_rho": np.nan,
        "tf_pval": np.nan,
        "tf_rho_conditional": np.nan,
        "tf_pval_conditional": np.nan,
        "sf_rho": np.nan,
        "sf_pval": np.nan,
        "delta_usage": np.nan,
        "qc_ok": False,
        "n_samples": 0,
    }
    if reg_tx not in expr_df_indexed.index or target_tx not in expr_df_indexed.index:
        return result
    if target_tx not in usage_df_indexed.index:
        return result

    reg_expr = np.array(expr_df_indexed.loc[reg_tx].values, dtype=np.float64)
    target_expr = np.array(expr_df_indexed.loc[target_tx].values, dtype=np.float64)
    target_usage = np.array(usage_df_indexed.loc[target_tx].values, dtype=np.float64)
    target_reliable = np.array(
        reliability_df_indexed.loc[target_tx].values, dtype=np.float64
    )

    if reg_expr.ndim > 1:
        reg_expr = reg_expr[0]
    if target_expr.ndim > 1:
        target_expr = target_expr[0]
    if target_usage.ndim > 1:
        target_usage = target_usage[0]
    if target_reliable.ndim > 1:
        target_reliable = target_reliable[0]

    valid = (
        (target_reliable > 0)
        & np.isfinite(reg_expr)
        & np.isfinite(target_expr)
        & np.isfinite(target_usage)
    )
    if valid.sum() < 20:
        return result

    result["n_samples"] = int(valid.sum())
    result["qc_ok"] = True

    tf_rho, tf_pval = spearmanr(reg_expr[valid], target_expr[valid])
    result["tf_rho"] = tf_rho
    result["tf_pval"] = tf_pval

    if sf_expr_matrix is not None and len(sf_expr_matrix.columns) > 0:
        valid_indices = np.where(valid)[0]
        sf_data = sf_expr_matrix.iloc[valid_indices].values
        sf_var = np.var(sf_data, axis=0)
        sf_cols_valid = sf_var > 1e-6
        if sf_cols_valid.sum() > 0:
            sf_data_filtered = sf_data[:, sf_cols_valid]
            lr = LinearRegression()
            lr.fit(sf_data_filtered, target_expr[valid])
            target_residuals = target_expr[valid] - lr.predict(sf_data_filtered)
            tf_rho_cond, tf_pval_cond = spearmanr(reg_expr[valid], target_residuals)
            result["tf_rho_conditional"] = tf_rho_cond
            result["tf_pval_conditional"] = tf_pval_cond
        else:
            result["tf_rho_conditional"] = tf_rho
            result["tf_pval_conditional"] = tf_pval
    else:
        result["tf_rho_conditional"] = tf_rho
        result["tf_pval_conditional"] = tf_pval

    sf_rho, sf_pval = spearmanr(reg_expr[valid], target_usage[valid])
    result["sf_rho"] = sf_rho
    result["sf_pval"] = sf_pval

    q25, q75 = np.percentile(reg_expr[valid], [25, 75])
    low_mask = reg_expr[valid] <= q25
    high_mask = reg_expr[valid] >= q75
    if low_mask.sum() > 0 and high_mask.sum() > 0:
        result["delta_usage"] = float(
            np.mean(target_usage[valid][high_mask])
            - np.mean(target_usage[valid][low_mask])
        )
    else:
        result["delta_usage"] = 0
    return result


def compute_sf_splicing_evidence(
    sf_tx: str,
    target_gene: str,
    target_tx_set,
    mean_importance_sum: float,
    mean_importance_max: float,
    median_importance_sum: float,
    median_importance_max: float,
    sf_gene: str,
    n_target_tx: int,
    sf_expr_dict: dict,
    usage_df_indexed: pd.DataFrame,
    reliability_df_indexed: pd.DataFrame,
    gene_to_all_transcripts: dict,
    transcript_mean_expr: dict,
    rho_min: float = 0.3,
    q_min: float = 0.05,
    du_min: float = 0.1,
    min_isoforms: int = 2,
    top_m: int = 3,
    eps: float = 1e-6,
) -> dict:
    """Compute SF splicing evidence for a SF transcript - target gene pair.

    Args:
        sf_tx: SF transcript ID.
        target_gene: Target gene ID.
        target_tx_set: Set or iterable of target transcript IDs from Network 3.
        mean_importance_sum: Sum of mean importances across target transcripts.
        mean_importance_max: Max mean importance across target transcripts.
        median_importance_sum: Sum of median importances across target transcripts.
        median_importance_max: Max median importance across target transcripts.
        sf_gene: SF gene ID.
        n_target_tx: Number of target transcripts in network.
        sf_expr_dict: Dict mapping SF transcript ID to expression array.
        usage_df_indexed: Isoform usage DataFrame indexed by transcript ID.
        reliability_df_indexed: Usage reliability DataFrame indexed by transcript ID.
        gene_to_all_transcripts: Dict mapping gene ID to list of transcript IDs.
        transcript_mean_expr: Dict mapping transcript ID to mean expression.
        rho_min: Minimum Spearman correlation threshold.
        q_min: Maximum FDR-adjusted p-value threshold.
        du_min: Minimum delta usage threshold.
        min_isoforms: Minimum number of isoforms required for splicing claims.
        top_m: Number of top expressed transcripts to test.
        eps: Small constant to avoid division by zero.

    Returns:
        Dictionary with SF splicing evidence metrics.
    """
    dom_min = 0.5

    result = {
        "sf_tx": sf_tx,
        "sf_gene": sf_gene,
        "target_gene": target_gene,
        "target_tx_set": (
            ",".join(sorted(target_tx_set))
            if isinstance(target_tx_set, set)
            else str(target_tx_set)
        ),
        "mean_importance_sum": mean_importance_sum,
        "mean_importance_max": mean_importance_max,
        "median_importance_sum": median_importance_sum,
        "median_importance_max": median_importance_max,
        "n_target_tx_in_net": n_target_tx,
        "n_isoforms_in_gene": 0,
        "n_transcripts_tested": 0,
        "n_sig": 0,
        "push_pull": False,
        "dominance": 0.0,
        "usage_reliable": False,
        "best_rho": 0.0,
        "best_tx": None,
        "sf_category": "sf_ambiguous",
    }

    if sf_tx not in sf_expr_dict:
        return result

    sf_vals = sf_expr_dict[sf_tx]
    all_gene_txs = gene_to_all_transcripts.get(target_gene, [])
    all_gene_txs_in_data = [
        tx for tx in all_gene_txs if tx in usage_df_indexed.index
    ]
    result["n_isoforms_in_gene"] = len(all_gene_txs_in_data)

    if len(all_gene_txs_in_data) < min_isoforms:
        return result

    txs_to_test = set(tx for tx in target_tx_set if tx in usage_df_indexed.index)
    gene_txs_sorted = sorted(
        all_gene_txs_in_data,
        key=lambda tx: transcript_mean_expr.get(tx, 0),
        reverse=True,
    )
    for tx in gene_txs_sorted[:top_m]:
        txs_to_test.add(tx)
    txs_to_test = list(txs_to_test)

    if len(txs_to_test) == 0:
        return result

    result["n_transcripts_tested"] = len(txs_to_test)
    result["usage_reliable"] = True

    correlations = []
    for tx in txs_to_test:
        tx_usage = np.array(usage_df_indexed.loc[tx].values, dtype=np.float64)
        tx_reliable = np.array(
            reliability_df_indexed.loc[tx].values, dtype=np.float64
        )
        if tx_usage.ndim > 1:
            tx_usage = tx_usage[0]
        if tx_reliable.ndim > 1:
            tx_reliable = tx_reliable[0]

        valid = (
            (tx_reliable > 0) & np.isfinite(sf_vals) & np.isfinite(tx_usage)
        )
        if valid.sum() < 10:
            continue

        rho, pval = spearmanr(sf_vals[valid], tx_usage[valid])
        if not np.isfinite(rho):
            continue

        q25, q75 = np.percentile(sf_vals[valid], [25, 75])
        low_mask = sf_vals[valid] <= q25
        high_mask = sf_vals[valid] >= q75
        delta_usage = 0.0
        if low_mask.sum() > 0 and high_mask.sum() > 0:
            delta_usage = float(
                np.mean(tx_usage[valid][high_mask])
                - np.mean(tx_usage[valid][low_mask])
            )

        correlations.append(
            {"tx": tx, "rho": rho, "pval": pval, "delta_usage": delta_usage}
        )

    if len(correlations) == 0:
        result["sf_category"] = "sf_expression_associated"
        return result

    _, qvals, _, _ = multipletests(
        [c["pval"] for c in correlations], method="fdr_bh"
    )
    for i, c in enumerate(correlations):
        c["qval"] = qvals[i]

    sig_corrs = [
        c
        for c in correlations
        if c["qval"] <= q_min
        and abs(c["rho"]) >= rho_min
        and abs(c["delta_usage"]) >= du_min
    ]
    n_sig = len(sig_corrs)
    result["n_sig"] = n_sig

    if n_sig >= 2:
        rhos = [c["rho"] for c in sig_corrs]
        result["push_pull"] = bool(max(rhos) > 0 and min(rhos) < 0)

    if n_sig > 0:
        abs_rhos = [abs(c["rho"]) for c in sig_corrs]
        result["dominance"] = max(abs_rhos) / (sum(abs_rhos) + eps)
        best_idx = int(np.argmax(abs_rhos))
        result["best_rho"] = sig_corrs[best_idx]["rho"]
        result["best_tx"] = sig_corrs[best_idx]["tx"]

    if n_sig >= 2:
        if result["push_pull"] or result["dominance"] >= dom_min:
            result["sf_category"] = "sf_splicing_supported_specific"
        else:
            result["sf_category"] = "sf_splicing_supported_broad"
    else:
        result["sf_category"] = "sf_expression_associated"

    return result


def _run_pipeline(
    canonical_grn: pd.DataFrame,
    as_source_grn: pd.DataFrame,
    fully_as_grn: pd.DataFrame,
    regulator_list: pd.DataFrame,
    transcript_tpm: pd.DataFrame,
    sample_cols: list,
    tx2gene: dict,
    min_frequency: int = 10,
    importance_percentile: float = 0.7,
    r_iso: float = 1.5,
    r_gene: float = 1.5,
    r_eq: float = 1.2,
    rho_min: float = 0.3,
    q_min: float = 0.05,
    du_min: float = 0.1,
    gene_tpm_min: float = 1.0,
    eps: float = 1e-6,
) -> dict:
    """Run the full edge categorization pipeline.

    Args:
        canonical_grn: Network 1 (canonical, gene-level) edge DataFrame.
        as_source_grn: Network 2 (AS-aware source) edge DataFrame.
        fully_as_grn: Network 3 (fully AS-aware) edge DataFrame.
        regulator_list: Regulator list with ``Regulator_type`` column.
        transcript_tpm: Transcript TPM DataFrame with ``transcript_id``
            and ``gene_id`` columns.
        sample_cols: List of sample column names in ``transcript_tpm``.
        tx2gene: Dict mapping transcript ID to gene ID.
        min_frequency: Minimum frequency threshold for filtering.
        importance_percentile: Importance percentile threshold for filtering.
        r_iso: Fold-change threshold for isoform-specific classification.
        r_gene: Fold-change threshold for gene-specific classification.
        r_eq: Equivalence band threshold.
        rho_min: Minimum Spearman correlation for evidence.
        q_min: Maximum q-value for evidence.
        du_min: Minimum delta usage for evidence.
        gene_tpm_min: Minimum gene TPM for reliability.
        eps: Small constant to avoid division by zero.

    Returns:
        Dictionary with keys: ``net1_filtered``, ``net2_filtered``,
        ``net3_filtered``, ``set_a``, ``set_b``, ``set_c``, ``set_d``,
        ``set_b_unpacked``, ``set_c_unpacked``.
    """
    net1_raw = canonical_grn.rename(
        columns={"source": "source_gene", "target": "target_gene"}
    )
    net2_raw = as_source_grn.rename(
        columns={"source": "source_transcript", "target": "target_gene"}
    )
    net2_raw["source_gene"] = net2_raw["source_transcript"].map(tx2gene)

    net3_raw = fully_as_grn.rename(
        columns={"source": "source_transcript", "target": "target_transcript"}
    )
    net3_raw["source_gene"] = net3_raw["source_transcript"].map(tx2gene)
    net3_raw["target_gene"] = net3_raw["target_transcript"].map(tx2gene)

    net1_filtered = filter_edges(net1_raw, min_frequency, importance_percentile)
    net2_filtered = filter_edges(net2_raw, min_frequency, importance_percentile)
    net3_filtered = filter_edges(net3_raw, min_frequency, importance_percentile)

    tx_to_regtype = dict(
        zip(regulator_list["Transcript stable ID"], regulator_list["Regulator_type"])
    )
    gene_to_regtype = (
        regulator_list.groupby("Gene stable ID")["Regulator_type"].first().to_dict()
    )

    net1_filtered = net1_filtered.copy()
    net1_filtered["reg_type"] = net1_filtered["source_gene"].map(gene_to_regtype)
    net2_filtered = net2_filtered.copy()
    net2_filtered["reg_type"] = net2_filtered["source_transcript"].map(tx_to_regtype)
    net3_filtered = net3_filtered.copy()
    net3_filtered["reg_type"] = net3_filtered["source_transcript"].map(tx_to_regtype)

    net3_tf = net3_filtered[net3_filtered["reg_type"] == "TF"].copy()
    net3_sf = net3_filtered[net3_filtered["reg_type"] == "SF"].copy()
    net3_tfsf = net3_filtered[net3_filtered["reg_type"] == "TF_SF"].copy()

    # Compute isoform usage and reliability
    tpm_numeric = transcript_tpm[sample_cols].astype(float)
    tpm_with_gene = tpm_numeric.copy()
    tpm_with_gene["_gene_id_"] = transcript_tpm["gene_id"].values
    gene_totals = tpm_with_gene.groupby("_gene_id_")[sample_cols].transform("sum")
    usage_values = tpm_numeric / (gene_totals + eps)

    usage_df = transcript_tpm[["transcript_id", "gene_id"]].copy()
    usage_df[sample_cols] = usage_values
    reliability_df = (gene_totals >= gene_tpm_min).astype(float)

    usage_df_indexed = usage_df.set_index("transcript_id")[sample_cols]
    reliability_df_indexed = reliability_df.copy()
    reliability_df_indexed.index = transcript_tpm["transcript_id"].values
    reliability_df_indexed = reliability_df_indexed[sample_cols]

    gene_to_all_transcripts = (
        usage_df.groupby("gene_id")["transcript_id"].apply(list).to_dict()
    )
    transcript_mean_expr = (
        transcript_tpm.set_index("transcript_id")[sample_cols].mean(axis=1).to_dict()
    )
    expr_df_indexed = transcript_tpm.set_index("transcript_id")[sample_cols]

    # Set A: source resolution (TF + TF_SF regulators)
    net1_tf = net1_filtered[net1_filtered["reg_type"].isin(["TF", "TF_SF"])].copy()
    net2_tf = net2_filtered[net2_filtered["reg_type"].isin(["TF", "TF_SF"])].copy()

    net1_tf["edge_gg"] = net1_tf["source_gene"] + "|" + net1_tf["target_gene"]
    net2_tf["edge_gg"] = net2_tf["source_gene"] + "|" + net2_tf["target_gene"]

    net1_agg = (
        net1_tf.groupby("edge_gg")
        .agg({"mean_importance": "max", "median_importance": "max"})
        .reset_index()
    )
    net1_mean_imp = dict(zip(net1_agg["edge_gg"], net1_agg["mean_importance"]))
    net1_median_imp = dict(zip(net1_agg["edge_gg"], net1_agg["median_importance"]))

    net2_agg = net2_tf.loc[net2_tf.groupby("edge_gg")["mean_importance"].idxmax()]
    net2_mean_imp = dict(zip(net2_agg["edge_gg"], net2_agg["mean_importance"]))
    net2_median_imp = dict(zip(net2_agg["edge_gg"], net2_agg["median_importance"]))
    net2_best_tx = dict(zip(net2_agg["edge_gg"], net2_agg["source_transcript"]))

    all_edges_t1 = set(net1_mean_imp.keys()) | set(net2_mean_imp.keys())
    set_a_rows = [
        categorize_source_resolution(
            e, net1_mean_imp, net2_mean_imp, net1_median_imp, net2_median_imp,
            net2_best_tx, r_iso, r_gene, r_eq, eps
        )
        for e in all_edges_t1
    ]
    set_a = pd.DataFrame(set_a_rows)
    set_a["max_median"] = set_a[["S1_median", "S2_median"]].max(axis=1)
    set_a = set_a.sort_values("max_median", ascending=False)

    # Set D: TF_SF joint categorization
    if len(net3_tfsf) > 0:
        sf_tx_list = regulator_list[
            regulator_list["Regulator_type"] == "SF"
        ]["Transcript stable ID"].tolist()
        sf_tx_in_expr = [tx for tx in sf_tx_list if tx in expr_df_indexed.index]
        if sf_tx_in_expr:
            sf_expr_matrix = expr_df_indexed.loc[sf_tx_in_expr].T
            sf_expr_matrix = sf_expr_matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        else:
            sf_expr_matrix = None

        tfsf_valid = net3_tfsf[
            net3_tfsf["source_transcript"].isin(expr_df_indexed.index)
            & net3_tfsf["target_transcript"].isin(expr_df_indexed.index)
        ].copy()

        set_d_rows = []
        for _, row in tfsf_valid.iterrows():
            evidence = compute_tfsf_evidence(
                row["source_transcript"],
                row["target_transcript"],
                expr_df_indexed,
                usage_df_indexed,
                reliability_df_indexed,
                sf_expr_matrix,
                sample_cols,
                rho_min,
                q_min,
                du_min,
            )
            set_d_rows.append(
                {
                    "reg_tx": row["source_transcript"],
                    "reg_gene": row["source_gene"],
                    "target_tx": row["target_transcript"],
                    "target_gene": row["target_gene"],
                    "mean_importance": row["mean_importance"],
                    "median_importance": row["median_importance"],
                    **evidence,
                }
            )

        set_d_full = pd.DataFrame(set_d_rows) if set_d_rows else pd.DataFrame()

        if len(set_d_full) > 0:
            valid_tf = set_d_full["tf_pval_conditional"].notna() & (
                set_d_full["tf_pval_conditional"] > 0
            )
            valid_sf = set_d_full["sf_pval"].notna() & (set_d_full["sf_pval"] > 0)
            set_d_full["tf_q"] = np.nan
            set_d_full["sf_q"] = np.nan
            if valid_tf.sum() > 0:
                _, tf_qvals, _, _ = multipletests(
                    set_d_full.loc[valid_tf, "tf_pval_conditional"].values,
                    method="fdr_bh",
                )
                set_d_full.loc[valid_tf, "tf_q"] = tf_qvals
            if valid_sf.sum() > 0:
                _, sf_qvals, _, _ = multipletests(
                    set_d_full.loc[valid_sf, "sf_pval"].values, method="fdr_bh"
                )
                set_d_full.loc[valid_sf, "sf_q"] = sf_qvals

            set_d_full["tf_evidence"] = (
                (set_d_full["tf_q"] <= q_min)
                & (set_d_full["tf_rho_conditional"].abs() >= rho_min)
            ).fillna(False)
            set_d_full["sf_evidence"] = (
                (set_d_full["sf_q"] <= q_min)
                & (set_d_full["sf_rho"].abs() >= rho_min)
                & (set_d_full["delta_usage"].abs() >= du_min)
            ).fillna(False)

            def _categorize_tfsf(row):
                if not row["qc_ok"]:
                    return "tfsf_ambiguous"
                tf_ev, sf_ev = row["tf_evidence"], row["sf_evidence"]
                if sf_ev and not tf_ev:
                    return "tfsf_sf_like"
                elif tf_ev and not sf_ev:
                    return "tfsf_tf_like"
                elif tf_ev and sf_ev:
                    return "tfsf_joint"
                return "tfsf_ambiguous"

            set_d_full["tfsf_category"] = set_d_full.apply(_categorize_tfsf, axis=1)
            set_d_full = set_d_full.sort_values("median_importance", ascending=False)
        else:
            if len(net3_tfsf) > 0:
                set_d_full = net3_tfsf[
                    ["source_transcript", "source_gene", "target_transcript",
                     "target_gene", "mean_importance", "median_importance"]
                ].copy()
                set_d_full.columns = [
                    "reg_tx", "reg_gene", "target_tx", "target_gene",
                    "mean_importance", "median_importance"
                ]
                set_d_full["tfsf_category"] = "tfsf_ambiguous"
    else:
        set_d_full = pd.DataFrame()

    tfsf_tf_like = (
        set_d_full[set_d_full["tfsf_category"] == "tfsf_tf_like"].copy()
        if len(set_d_full) > 0
        else pd.DataFrame()
    )
    tfsf_sf_like = (
        set_d_full[set_d_full["tfsf_category"] == "tfsf_sf_like"].copy()
        if len(set_d_full) > 0
        else pd.DataFrame()
    )
    set_d = (
        set_d_full[set_d_full["tfsf_category"].isin(["tfsf_joint", "tfsf_ambiguous"])].copy()
        if len(set_d_full) > 0
        else pd.DataFrame()
    )

    # Set B: target resolution (TF edges)
    net2_for_t2 = net2_filtered[
        net2_filtered["reg_type"].isin(["TF", "TF_SF"])
    ].copy()
    net3_tf_only = net3_filtered[net3_filtered["reg_type"] == "TF"].copy()

    net2_for_t2["edge_tg"] = (
        net2_for_t2["source_transcript"] + "|" + net2_for_t2["target_gene"]
    )
    net2_tg_mean_imp = dict(
        zip(net2_for_t2["edge_tg"], net2_for_t2["mean_importance"])
    )
    net2_tg_median_imp = dict(
        zip(net2_for_t2["edge_tg"], net2_for_t2["median_importance"])
    )

    net3_tf_only["edge_tg"] = (
        net3_tf_only["source_transcript"] + "|" + net3_tf_only["target_gene"]
    )
    net3_tg_agg = (
        net3_tf_only.groupby("edge_tg")
        .agg(
            {
                "mean_importance": ["sum", "max"],
                "median_importance": ["sum", "max"],
                "source_gene": "first",
                "target_gene": "first",
                "target_transcript": lambda x: ",".join(sorted(set(x))),
            }
        )
        .reset_index()
    )
    net3_tg_agg.columns = [
        "edge_tg", "S3_mean_sum", "S3_mean_max",
        "S3_median_sum", "S3_median_max",
        "source_gene", "target_gene", "target_tx_set",
    ]
    net3_tg_agg["dominance"] = net3_tg_agg["S3_mean_max"] / (
        net3_tg_agg["S3_mean_sum"] + eps
    )

    net3_tg_mean_imp = dict(zip(net3_tg_agg["edge_tg"], net3_tg_agg["S3_mean_sum"]))
    net3_tg_mean_max = dict(zip(net3_tg_agg["edge_tg"], net3_tg_agg["S3_mean_max"]))
    net3_tg_median_imp = dict(
        zip(net3_tg_agg["edge_tg"], net3_tg_agg["S3_median_sum"])
    )
    net3_tg_dom = dict(zip(net3_tg_agg["edge_tg"], net3_tg_agg["dominance"]))

    all_edges_t2_tf = set(net2_tg_mean_imp.keys()) | set(net3_tg_mean_imp.keys())
    set_b_tf_rows = [
        categorize_target_resolution(
            e, net2_tg_mean_imp, net3_tg_mean_imp, net3_tg_mean_max, net3_tg_dom,
            net2_tg_median_imp, net3_tg_median_imp, tx2gene, r_iso, r_gene, r_eq, eps
        )
        for e in all_edges_t2_tf
    ]
    set_b_tf = pd.DataFrame(set_b_tf_rows)
    set_b_tf["reg_type"] = "TF"
    _net3_tx_map = dict(zip(net3_tg_agg["edge_tg"], net3_tg_agg["target_tx_set"]))
    set_b_tf["target_tx_set"] = set_b_tf.apply(
        lambda r: _net3_tx_map.get(r["regulator_tx"] + "|" + r["target_gene"], ""),
        axis=1,
    )

    if len(tfsf_tf_like) > 0:
        net2_tfsf = net2_filtered[net2_filtered["reg_type"] == "TF_SF"].copy()
        net2_tfsf["edge_tg"] = (
            net2_tfsf["source_transcript"] + "|" + net2_tfsf["target_gene"]
        )
        net2_tfsf_mean_imp = dict(
            zip(net2_tfsf["edge_tg"], net2_tfsf["mean_importance"])
        )
        net2_tfsf_median_imp = dict(
            zip(net2_tfsf["edge_tg"], net2_tfsf["median_importance"])
        )

        net3_tfsf_agg = (
            net3_tfsf.groupby(["source_transcript", "target_gene"])
            .agg(
                {
                    "mean_importance": ["sum", "max"],
                    "median_importance": ["sum", "max"],
                    "source_gene": "first",
                    "target_transcript": lambda x: ",".join(sorted(set(x))),
                }
            )
            .reset_index()
        )
        net3_tfsf_agg.columns = [
            "source_transcript", "target_gene",
            "S3_mean_sum", "S3_mean_max",
            "S3_median_sum", "S3_median_max",
            "source_gene", "target_tx_set",
        ]
        net3_tfsf_agg["edge_tg"] = (
            net3_tfsf_agg["source_transcript"] + "|" + net3_tfsf_agg["target_gene"]
        )
        net3_tfsf_agg["dominance"] = net3_tfsf_agg["S3_mean_max"] / (
            net3_tfsf_agg["S3_mean_sum"] + eps
        )

        tfsf_tg_mean_imp = dict(
            zip(net3_tfsf_agg["edge_tg"], net3_tfsf_agg["S3_mean_sum"])
        )
        tfsf_tg_mean_max = dict(
            zip(net3_tfsf_agg["edge_tg"], net3_tfsf_agg["S3_mean_max"])
        )
        tfsf_tg_median_imp = dict(
            zip(net3_tfsf_agg["edge_tg"], net3_tfsf_agg["S3_median_sum"])
        )
        tfsf_tg_dom = dict(zip(net3_tfsf_agg["edge_tg"], net3_tfsf_agg["dominance"]))

        tfsf_tf_like["edge_tg"] = (
            tfsf_tf_like["reg_tx"] + "|" + tfsf_tf_like["target_gene"]
        )
        tfsf_tf_like_edges = set(tfsf_tf_like["edge_tg"].unique())

        set_b_tfsf_rows = [
            categorize_target_resolution(
                e, net2_tfsf_mean_imp, tfsf_tg_mean_imp, tfsf_tg_mean_max, tfsf_tg_dom,
                net2_tfsf_median_imp, tfsf_tg_median_imp, tx2gene, r_iso, r_gene, r_eq, eps
            )
            for e in tfsf_tf_like_edges
        ]
        set_b_tfsf = pd.DataFrame(set_b_tfsf_rows)
        set_b_tfsf["reg_type"] = "TF_SF"
        _tfsf_tx_map = dict(
            zip(
                net3_tfsf_agg["source_transcript"] + "|" + net3_tfsf_agg["target_gene"],
                net3_tfsf_agg["target_tx_set"],
            )
        )
        set_b_tfsf["target_tx_set"] = set_b_tfsf.apply(
            lambda r: _tfsf_tx_map.get(
                r["regulator_tx"] + "|" + r["target_gene"], ""
            ),
            axis=1,
        )
    else:
        set_b_tfsf = pd.DataFrame()

    set_b = pd.concat([set_b_tf, set_b_tfsf], ignore_index=True)
    set_b["max_median"] = set_b[["S2_median", "S3_median"]].max(axis=1)
    set_b = set_b.sort_values("max_median", ascending=False)

    # Set B Unpacked
    net3_tf_unpack = net3_tf_only[
        ["source_transcript", "source_gene", "target_transcript",
         "target_gene", "mean_importance", "median_importance", "frequency"]
    ].copy()
    net3_tf_unpack.columns = [
        "regulator_tx", "regulator_gene", "target_tx", "target_gene",
        "net3_mean_importance", "net3_median_importance", "net3_frequency",
    ]
    net3_tf_unpack["reg_type"] = "TF"

    if len(tfsf_tf_like) > 0:
        tfsf_tf_like_regtx_set = set(tfsf_tf_like["reg_tx"].unique())
        net3_tfsf_unpack = net3_tfsf[
            net3_tfsf["source_transcript"].isin(tfsf_tf_like_regtx_set)
        ][
            ["source_transcript", "source_gene", "target_transcript",
             "target_gene", "mean_importance", "median_importance", "frequency"]
        ].copy()
        net3_tfsf_unpack.columns = [
            "regulator_tx", "regulator_gene", "target_tx", "target_gene",
            "net3_mean_importance", "net3_median_importance", "net3_frequency",
        ]
        net3_tfsf_unpack["reg_type"] = "TF_SF"
        set_b_unpacked = pd.concat(
            [net3_tf_unpack, net3_tfsf_unpack], ignore_index=True
        )
    else:
        set_b_unpacked = net3_tf_unpack.copy()

    set_b_unpacked["edge_tg"] = (
        set_b_unpacked["regulator_tx"] + "|" + set_b_unpacked["target_gene"]
    )
    set_b_key = set_b[["regulator_tx", "target_gene", "target_category"]].copy()
    set_b_key["edge_tg"] = set_b_key["regulator_tx"] + "|" + set_b_key["target_gene"]
    set_b_unpacked = set_b_unpacked.merge(
        set_b_key[["edge_tg", "target_category"]], on="edge_tg", how="inner"
    )
    set_b_unpacked["target_rank_within_edge"] = (
        set_b_unpacked.groupby("edge_tg")["net3_mean_importance"]
        .rank(ascending=False, method="min")
        .astype(int)
    )
    set_b_unpacked = set_b_unpacked.sort_values(
        "net3_median_importance", ascending=False
    )
    set_b_unpacked = set_b_unpacked.drop(columns=["edge_tg"])

    # Set C: SF splicing evidence
    sf_edges = (
        net3_sf.groupby(["source_transcript", "source_gene", "target_gene"])
        .agg(
            {
                "mean_importance": ["sum", "max", "count"],
                "median_importance": ["sum", "max"],
                "target_transcript": lambda x: set(x),
            }
        )
        .reset_index()
    )
    sf_edges.columns = [
        "sf_tx", "sf_gene", "target_gene",
        "mean_importance_sum", "mean_importance_max", "n_target_tx",
        "median_importance_sum", "median_importance_max", "target_tx_set",
    ]

    sf_transcripts = sf_edges["sf_tx"].unique()
    sf_expr_data = transcript_tpm[
        transcript_tpm["transcript_id"].isin(sf_transcripts)
    ].copy()
    sf_expr_data = sf_expr_data.set_index("transcript_id")[sample_cols]
    sf_expr_dict = {
        tx: np.array(sf_expr_data.loc[tx].values, dtype=np.float64)
        for tx in sf_expr_data.index
    }

    set_c_sf_rows = []
    for _, row in sf_edges.iterrows():
        result = compute_sf_splicing_evidence(
            row["sf_tx"],
            row["target_gene"],
            row["target_tx_set"],
            row["mean_importance_sum"],
            row["mean_importance_max"],
            row["median_importance_sum"],
            row["median_importance_max"],
            row["sf_gene"],
            row["n_target_tx"],
            sf_expr_dict,
            usage_df_indexed,
            reliability_df_indexed,
            gene_to_all_transcripts,
            transcript_mean_expr,
            rho_min,
            q_min,
            du_min,
        )
        set_c_sf_rows.append(result)

    set_c_sf_base = pd.DataFrame(set_c_sf_rows) if set_c_sf_rows else pd.DataFrame()

    if len(tfsf_sf_like) > 0 and len(set_c_sf_base) > 0:
        tfsf_sf_rows = []
        for _, row in tfsf_sf_like.iterrows():
            r = {
                "sf_tx": row["reg_tx"],
                "sf_gene": row["reg_gene"],
                "target_gene": row["target_gene"],
                "target_tx_set": row.get("target_tx", ""),
                "mean_importance_sum": row.get("mean_importance", 0),
                "mean_importance_max": row.get("mean_importance", 0),
                "median_importance_sum": row.get("median_importance", 0),
                "median_importance_max": row.get("median_importance", 0),
                "sf_category": row.get("tfsf_category", "tfsf_sf_like"),
                "reg_type": "TF_SF",
            }
            tfsf_sf_rows.append(r)
        tfsf_sf_df = pd.DataFrame(tfsf_sf_rows)
        set_c = pd.concat([set_c_sf_base, tfsf_sf_df], ignore_index=True)
    else:
        set_c = set_c_sf_base.copy()

    if len(set_c) > 0 and "reg_type" not in set_c.columns:
        set_c["reg_type"] = "SF"

    # Set C Unpacked
    set_c_unpacked_rows = []
    for _, row in net3_sf.iterrows():
        set_c_unpacked_rows.append(
            {
                "sf_tx": row["source_transcript"],
                "sf_gene": row["source_gene"],
                "target_tx": row["target_transcript"],
                "target_gene": row["target_gene"],
                "net3_mean_importance": row["mean_importance"],
                "net3_median_importance": row["median_importance"],
                "net3_frequency": row["frequency"],
                "reg_type": "SF",
            }
        )
    set_c_unpacked = pd.DataFrame(set_c_unpacked_rows) if set_c_unpacked_rows else pd.DataFrame()

    if len(set_c_unpacked) > 0 and len(set_c) > 0:
        set_c_unpacked["edge_sg"] = (
            set_c_unpacked["sf_tx"] + "|" + set_c_unpacked["target_gene"]
        )
        set_c_key = set_c[["sf_tx", "target_gene", "sf_category"]].copy()
        set_c_key["edge_sg"] = set_c_key["sf_tx"] + "|" + set_c_key["target_gene"]
        set_c_unpacked = set_c_unpacked.merge(
            set_c_key[["edge_sg", "sf_category"]], on="edge_sg", how="inner"
        )
        set_c_unpacked["target_rank_within_edge"] = (
            set_c_unpacked.groupby("edge_sg")["net3_mean_importance"]
            .rank(ascending=False, method="min")
            .astype(int)
        )
        set_c_unpacked = set_c_unpacked.drop(columns=["edge_sg"])

    return {
        "net1_filtered": net1_filtered,
        "net2_filtered": net2_filtered,
        "net3_filtered": net3_filtered,
        "set_a": set_a,
        "set_b": set_b,
        "set_c": set_c,
        "set_d": set_d,
        "set_b_unpacked": set_b_unpacked,
        "set_c_unpacked": set_c_unpacked,
    }
