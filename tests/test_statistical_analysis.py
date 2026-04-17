"""Tests for alternet.statistical_analysis module."""

import numpy as np
import pandas as pd
import pytest

from alternet.statistical_analysis import (
    compute_section_A,
    compute_section_B,
    compute_section_C,
)


@pytest.fixture
def sample_set_b():
    """Minimal Set B DataFrame for section A tests."""
    return pd.DataFrame(
        {
            "regulator_tx": ["ENST001", "ENST001", "ENST001"],
            "regulator_gene": ["ENSG001", "ENSG001", "ENSG001"],
            "target_gene": ["ENSG002", "ENSG002", "ENSG003"],
            "target_category": [
                "target_isoform_specific",
                "target_gene_specific",
                "target_isoform_specific",
            ],
            "tgt_n_isoforms": [3, 2, 4],
            "tgt_dominance": [0.6, 0.7, 0.5],
            "target_tx_set": ["ENST003,ENST004", "", "ENST005"],
        }
    )


@pytest.fixture
def sample_set_c():
    """Minimal Set C DataFrame for section B tests."""
    return pd.DataFrame(
        {
            "sf_tx": ["ENST004", "ENST004"],
            "sf_gene": ["ENSG003", "ENSG003"],
            "target_gene": ["ENSG001", "ENSG002"],
            "sf_category": [
                "sf_splicing_supported_specific",
                "sf_expression_associated",
            ],
            "tgt_n_isoforms": [2, 3],
            "dominance": [0.7, 0.5],
            "target_tx_set": ["ENST001,ENST002", "ENST003"],
        }
    )


@pytest.fixture
def sample_set_a():
    """Minimal Set A DataFrame for section C tests."""
    return pd.DataFrame(
        {
            "regulator_gene": ["ENSG001", "ENSG001"],
            "target_gene": ["ENSG002", "ENSG003"],
            "source_category": ["source_isoform_specific", "source_gene_specific"],
            "best_tx": ["ENST001", "ENST002"],
            "reg_appris": ["PRINCIPAL:1", "ALTERNATIVE:1"],
        }
    )


def test_compute_section_A(sample_set_b):
    data = {"t2": sample_set_b}
    result = compute_section_A(data)
    assert "stats" in result
    assert "plot_data" in result
    assert "n_edges" in result["stats"]
    assert result["stats"]["n_edges"] == 3


def test_compute_section_B(sample_set_c):
    data = {"t3": sample_set_c}
    result = compute_section_B(data)
    assert "stats" in result
    assert result["stats"]["n_supported"] == 1


def test_compute_section_C(sample_set_a):
    tx_appris_label = {
        "ENST001": "PRINCIPAL:1",
        "ENST002": "ALTERNATIVE:1",
    }
    tx_trifid = {"ENST001": 0.9, "ENST002": 0.4}
    gene2txs = {"ENSG001": {"ENST001", "ENST002"}}
    data = {"t1": sample_set_a}
    result = compute_section_C(data, tx_appris_label, tx_trifid, gene2txs)
    assert "stats" in result
    assert "n_edges" in result["stats"]
