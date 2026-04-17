"""Tests for alternet.edge_categorization module."""

import numpy as np
import pandas as pd
import pytest

from alternet.edge_categorization import (
    categorize_source_resolution,
    categorize_target_resolution,
    compute_sf_splicing_evidence,
    compute_tfsf_evidence,
    filter_edges,
)


def test_filter_edges():
    df = pd.DataFrame(
        {
            "frequency": [5, 15, 20, 8, 25],
            "median_importance": [0.1, 0.8, 0.9, 0.3, 0.95],
        }
    )
    result = filter_edges(df, min_frequency=10, importance_percentile=0.5)
    assert (result["frequency"] >= 10).all()
    assert len(result) < len(df)


def test_categorize_source_resolution_isoform_specific():
    result = categorize_source_resolution(
        "ENSG001|ENSG002",
        net1_mean={"ENSG001|ENSG002": 0.2},
        net2_mean={"ENSG001|ENSG002": 0.6},
        net1_median={"ENSG001|ENSG002": 0.15},
        net2_median={"ENSG001|ENSG002": 0.55},
        net2_best_tx={"ENSG001|ENSG002": "ENST001"},
        r_iso=1.5,
        r_gene=1.5,
        r_eq=1.2,
    )
    assert result["source_category"] == "source_isoform_specific"


def test_categorize_source_resolution_gene_specific():
    result = categorize_source_resolution(
        "ENSG001|ENSG002",
        net1_mean={"ENSG001|ENSG002": 0.6},
        net2_mean={"ENSG001|ENSG002": 0.2},
        net1_median={},
        net2_median={},
        net2_best_tx={},
        r_iso=1.5,
        r_gene=1.5,
        r_eq=1.2,
    )
    assert result["source_category"] == "source_gene_specific"


def test_categorize_target_resolution_isoform_specific():
    tx2gene = {"ENST001": "ENSG001"}
    result = categorize_target_resolution(
        "ENST001|ENSG002",
        net2_mean={"ENST001|ENSG002": 0.2},
        net3_mean={"ENST001|ENSG002": 0.6},
        net3_max={"ENST001|ENSG002": 0.6},
        net3_dom={"ENST001|ENSG002": 1.0},
        net2_median={"ENST001|ENSG002": 0.15},
        net3_median={"ENST001|ENSG002": 0.55},
        tx2gene=tx2gene,
    )
    assert result["target_category"] == "target_isoform_specific"


def test_compute_tfsf_evidence_missing_tx():
    expr = pd.DataFrame({"S1": [1.0, 2.0]}, index=["ENST001", "ENST002"])
    usage = pd.DataFrame({"S1": [0.6, 0.4]}, index=["ENST001", "ENST002"])
    reliability = pd.DataFrame({"S1": [1.0, 1.0]}, index=["ENST001", "ENST002"])
    result = compute_tfsf_evidence(
        "MISSING", "ENST002", expr, usage, reliability, None, ["S1"]
    )
    assert not result["qc_ok"]


def test_compute_sf_splicing_evidence_insufficient_isoforms():
    usage = pd.DataFrame(
        {"S1": [0.7], "S2": [0.3]}, index=["ENST001"]
    )
    reliability = pd.DataFrame(
        {"S1": [1.0], "S2": [1.0]}, index=["ENST001"]
    )
    result = compute_sf_splicing_evidence(
        sf_tx="ENST_SF",
        target_gene="ENSG001",
        target_tx_set={"ENST001"},
        mean_importance_sum=1.0,
        mean_importance_max=1.0,
        median_importance_sum=0.9,
        median_importance_max=0.9,
        sf_gene="ENSG_SF",
        n_target_tx=1,
        sf_expr_dict={"ENST_SF": np.array([1.0, 2.0])},
        usage_df_indexed=usage,
        reliability_df_indexed=reliability,
        gene_to_all_transcripts={"ENSG001": ["ENST001"]},
        transcript_mean_expr={"ENST001": 5.0},
        min_isoforms=2,
    )
    assert result["sf_category"] == "sf_ambiguous"
