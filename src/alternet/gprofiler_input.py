"""g:Profiler input generation pipeline for AlterNet 2.0."""

import pandas as pd


def score_targets(
    edges: pd.DataFrame,
    target_col: str,
    importance_col: str = "median_importance",
) -> pd.Series:
    """Compute target scores as weighted in-degree.

    score(target) = sum of importance for all edges targeting it.

    Args:
        edges: Edge DataFrame.
        target_col: Column name containing target IDs.
        importance_col: Column name containing edge importance scores.

    Returns:
        Series mapping target ID to score, sorted descending.
    """
    scores = edges.groupby(target_col)[importance_col].sum()
    return scores.sort_values(ascending=False)


def project_tx_to_gene(
    tx_scores: pd.Series,
    tx2gene: dict,
    method: str = "max",
) -> tuple:
    """Project transcript-level scores to gene-level.

    Args:
        tx_scores: Series mapping transcript_id to score.
        tx2gene: Dict mapping transcript_id to gene_id.
        method: Aggregation method: ``'max'`` (default) or ``'sum'``.

    Returns:
        Tuple of (gene_scores Series, rep_tx dict mapping gene_id to
        representative transcript).
    """
    df = pd.DataFrame(
        {"transcript_id": tx_scores.index, "score": tx_scores.values}
    )
    df["gene_id"] = df["transcript_id"].map(tx2gene)
    df = df.dropna(subset=["gene_id"])

    if method == "max":
        idx = df.groupby("gene_id")["score"].idxmax()
        result = df.loc[idx].set_index("gene_id")
        gene_scores = result["score"].sort_values(ascending=False)
        rep_tx = result["transcript_id"].to_dict()
    elif method == "sum":
        gene_scores = df.groupby("gene_id")["score"].sum().sort_values(ascending=False)
        rep_tx = {}
    else:
        raise ValueError(f"Unknown method: {method}")

    return gene_scores, rep_tx


def build_target_list(
    edges: pd.DataFrame,
    target_col: str,
    importance_col: str = "median_importance",
    target_type: str = "gene",
    tx2gene: dict = None,
    gene2symbol: dict = None,
    K: int = 500,
    projection_method: str = "max",
) -> tuple:
    """Build a ranked target list from an edge set for g:Profiler input.

    Args:
        edges: Edge DataFrame.
        target_col: Column name containing target IDs.
        importance_col: Column name containing importance scores.
        target_type: ``'gene'`` or ``'transcript'``.
        tx2gene: Transcript to gene mapping (required if target_type='transcript').
        gene2symbol: Gene ID to symbol mapping.
        K: Number of top genes to return.
        projection_method: ``'max'`` or ``'sum'`` for transcript projection.

    Returns:
        Tuple of (top_df DataFrame, symbols_list, meta_dict).
    """
    if len(edges) == 0:
        return pd.DataFrame(), [], {"n_edges": 0}

    target_scores = score_targets(edges, target_col, importance_col)

    if target_type == "transcript":
        if tx2gene is None:
            raise ValueError("tx2gene mapping required for transcript targets")
        gene_scores, rep_tx = project_tx_to_gene(
            target_scores, tx2gene, method=projection_method
        )
    else:
        gene_scores = target_scores
        rep_tx = {}

    top_genes = gene_scores.head(K)

    top_df = pd.DataFrame(
        {
            "rank": range(1, len(top_genes) + 1),
            "gene_id": top_genes.index,
            "score": top_genes.values,
        }
    )

    if gene2symbol is not None:
        top_df["gene_symbol"] = top_df["gene_id"].map(gene2symbol)
        top_df = top_df.dropna(subset=["gene_symbol"])
        gene_symbols = top_df["gene_symbol"].tolist()
    else:
        gene_symbols = top_df["gene_id"].tolist()

    if len(rep_tx) > 0:
        top_df["rep_transcript"] = top_df["gene_id"].map(rep_tx)

    metadata = {
        "n_edges": len(edges),
        "n_targets": len(target_scores),
        "n_genes": len(gene_scores),
        "n_top_symbols": len(gene_symbols),
        "projection_method": projection_method if target_type == "transcript" else "none",
    }

    return top_df, gene_symbols, metadata


def _run_pipeline(
    set_a: pd.DataFrame,
    set_b: pd.DataFrame,
    net1: pd.DataFrame,
    net2: pd.DataFrame,
    net3: pd.DataFrame,
    tx2gene: dict,
    gene2symbol: dict,
    top_k: int = 500,
) -> dict:
    """Build L1, L2, L3 gene lists for g:Profiler enrichment analysis.

    Args:
        set_a: Set A (source resolution) DataFrame from plausibility filtering.
        set_b: Set B (target resolution) DataFrame from plausibility filtering.
        net1: Network 1 filtered edge DataFrame.
        net2: Network 2 filtered edge DataFrame.
        net3: Network 3 filtered edge DataFrame.
        tx2gene: Dict mapping transcript ID to gene ID.
        gene2symbol: Dict mapping gene ID to gene symbol.
        top_k: Number of top genes to include in each list.

    Returns:
        Dictionary with keys: ``l1_top``, ``l1_symbols``, ``l2_top``,
        ``l2_symbols``, ``l3_top``, ``l3_symbols``.
    """
    l1_top, l1_symbols, _ = build_target_list(
        edges=net1,
        target_col="target_gene",
        importance_col="mean_importance",
        target_type="gene",
        gene2symbol=gene2symbol,
        K=top_k,
    )

    l2_edges = set_a[set_a["source_category"] == "source_isoform_specific"].copy()
    importance_col_l2 = "S2_median" if "S2_median" in l2_edges.columns else "median_importance"
    l2_top, l2_symbols, _ = build_target_list(
        edges=l2_edges,
        target_col="target_gene",
        importance_col=importance_col_l2,
        target_type="gene",
        gene2symbol=gene2symbol,
        K=top_k,
    )

    l3_table2 = set_b[set_b["target_category"] == "target_isoform_specific"].copy()
    if "reg_type" in net3.columns:
        net3_tf_like = net3[net3["reg_type"].isin(["TF", "TF_SF"])].copy()
    else:
        net3_tf_like = net3.copy()

    if "target_gene" not in net3_tf_like.columns and "target_transcript" in net3_tf_like.columns:
        net3_tf_like["target_gene"] = net3_tf_like["target_transcript"].map(tx2gene)

    target_col_net3 = (
        "target_transcript" if "target_transcript" in net3_tf_like.columns else "target"
    )
    source_col_net3 = (
        "source_transcript" if "source_transcript" in net3_tf_like.columns else "source"
    )

    if len(l3_table2) > 0 and len(net3_tf_like) > 0:
        reg_col = (
            "regulator_tx" if "regulator_tx" in l3_table2.columns else "source_transcript"
        )
        if reg_col in l3_table2.columns:
            l3_pairs = set(zip(l3_table2[reg_col], l3_table2["target_gene"]))
            net3_tf_like["pair"] = list(
                zip(net3_tf_like[source_col_net3], net3_tf_like["target_gene"])
            )
            l3_edges = net3_tf_like[net3_tf_like["pair"].isin(l3_pairs)].copy()
        else:
            l3_edges = net3_tf_like.copy()
    else:
        l3_edges = net3_tf_like.copy()

    l3_top, l3_symbols, _ = build_target_list(
        edges=l3_edges,
        target_col=target_col_net3,
        importance_col="mean_importance",
        target_type="transcript",
        tx2gene=tx2gene,
        gene2symbol=gene2symbol,
        K=top_k,
        projection_method="max",
    )

    return {
        "l1_top": l1_top,
        "l1_symbols": l1_symbols,
        "l2_top": l2_top,
        "l2_symbols": l2_symbols,
        "l3_top": l3_top,
        "l3_symbols": l3_symbols,
    }
