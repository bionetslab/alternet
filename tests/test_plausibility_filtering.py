"""Tests for alternet.plausibility_filtering module."""

import numpy as np
import pandas as pd
import pytest

from alternet.plausibility_filtering import (
    annotate_digger_raw,
    compute_dominance_metrics,
    filter_set_a,
    filter_set_b,
    filter_set_c,
    filter_set_d,
    get_principal_tx,
    get_shared_features,
    get_unique_features,
)


def test_compute_dominance_metrics(transcript_tpm_df, sample_cols):
    gene_dominance, gene_n_isoforms, tx_expr_share = compute_dominance_metrics(
        transcript_tpm_df, sample_cols
    )
    assert "ENSG001" in gene_dominance
    assert gene_dominance["ENSG001"] <= 1.0
    assert gene_n_isoforms["ENSG001"] == 2


def test_filter_set_a():
    set_a = pd.DataFrame(
        {
            "regulator_gene": ["ENSG001", "ENSG001"],
            "target_gene": ["ENSG002", "ENSG003"],
            "source_category": ["source_isoform_specific", "source_equivalent"],
            "best_tx": ["ENST001", "ENST001"],
            "S1_mean": [0.3, 0.5],
            "S2_mean": [0.9, 0.55],
            "S1_median": [0.25, 0.45],
            "S2_median": [0.85, 0.5],
        }
    )
    gene_dominance = {"ENSG001": 0.5}
    gene_n_isoforms = {"ENSG001": 2}
    result = filter_set_a(set_a, gene_dominance, gene_n_isoforms)
    assert "is_plausible" in result.columns
    assert len(result) == 2


def test_filter_set_b():
    set_b = pd.DataFrame(
        {
            "regulator_tx": ["ENST001"],
            "regulator_gene": ["ENSG001"],
            "target_gene": ["ENSG002"],
            "target_category": ["target_isoform_specific"],
            "S2_mean": [0.4],
            "S3_mean_sum": [1.0],
            "S2_median": [0.35],
            "S3_median": [0.9],
        }
    )
    gene_dominance = {"ENSG002": 0.6}
    gene_n_isoforms = {"ENSG002": 3}
    result = filter_set_b(set_b, gene_dominance, gene_n_isoforms)
    assert "is_plausible" in result.columns


def test_filter_set_c():
    set_c = pd.DataFrame(
        {
            "sf_tx": ["ENST004"],
            "sf_gene": ["ENSG003"],
            "target_gene": ["ENSG001"],
            "sf_category": ["sf_splicing_supported_specific"],
            "n_sig": [1],
            "push_pull": [False],
            "usage_reliable": [True],
        }
    )
    gene_n_isoforms = {"ENSG001": 1}
    result = filter_set_c(set_c, gene_n_isoforms, min_tx_per_gene=2)
    assert "is_plausible" in result.columns
    assert result["sf_category"].values[0] == "sf_expression_associated"


def test_filter_set_d():
    set_d = pd.DataFrame(
        {
            "reg_tx": ["ENST001"],
            "reg_gene": ["ENSG001"],
            "target_tx": ["ENST003"],
            "target_gene": ["ENSG002"],
            "tfsf_category": ["tfsf_joint"],
            "qc_ok": [True],
            "mean_importance": [0.8],
            "median_importance": [0.75],
        }
    )
    gene_dominance = {"ENSG002": 0.5}
    gene_n_isoforms = {"ENSG002": 2}
    result = filter_set_d(set_d, gene_dominance, gene_n_isoforms)
    assert "is_plausible" in result.columns


def test_get_principal_tx():
    digger_gene2tx = {"ENSG001": {"ENST001", "ENST002"}}
    tx_appris_label = {"ENST001": "PRINCIPAL:1", "ENST002": "ALTERNATIVE:1"}
    result = get_principal_tx("ENSG001", digger_gene2tx, tx_appris_label)
    assert result == "ENST001"


def test_get_unique_features():
    tx_exon_set = {"ENST001": {1, 2, 3}, "ENST002": {1, 2}}
    tx_domain_set = {"ENST001": {"PF001", "PF002"}, "ENST002": {"PF001"}}
    digger_gene2tx = {"ENSG001": {"ENST001", "ENST002"}}
    tx_appris_label = {"ENST001": "PRINCIPAL:1", "ENST002": "ALTERNATIVE:1"}
    ue, ud = get_unique_features(
        "ENST002", "ENSG001", tx_exon_set, tx_domain_set, digger_gene2tx, tx_appris_label
    )
    # ENST002 compared against principal ENST001: unique exons = {} (subset)
    assert isinstance(ue, str)
    assert isinstance(ud, str)


def test_get_shared_features():
    tx_exon_set = {"ENST001": {1, 2, 3}, "ENST002": {1, 2}}
    tx_domain_set = {"ENST001": {"PF001"}, "ENST002": {"PF001"}}
    digger_gene2tx = {"ENSG001": {"ENST001", "ENST002"}}
    se, sd = get_shared_features("ENSG001", digger_gene2tx, tx_exon_set, tx_domain_set)
    assert "1" in se
    assert "2" in se
    assert "3" not in se


def test_annotate_digger_raw():
    tx_exon_count = {"ENST001": 5}
    tx_domain_count = {"ENST001": 2}
    tx_domains = {"ENST001": "PF001,PF002"}
    ec, dc, ds = annotate_digger_raw("ENST001", tx_exon_count, tx_domain_count, tx_domains)
    assert ec == 5
    assert dc == 2
    assert ds == "PF001,PF002"
