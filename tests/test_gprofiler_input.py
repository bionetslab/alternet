"""Tests for alternet.gprofiler_input module."""

import pandas as pd
import pytest

from alternet.gprofiler_input import (
    build_target_list,
    project_tx_to_gene,
    score_targets,
)


def test_score_targets():
    edges = pd.DataFrame(
        {
            "target": ["ENSG001", "ENSG001", "ENSG002"],
            "median_importance": [0.5, 0.3, 0.8],
        }
    )
    result = score_targets(edges, "target", "median_importance")
    assert result["ENSG001"] == pytest.approx(0.8)
    assert result["ENSG002"] == pytest.approx(0.8)
    assert result.index[0] == "ENSG001" or result.index[0] == "ENSG002"


def test_project_tx_to_gene():
    tx_scores = pd.Series({"ENST001": 0.9, "ENST002": 0.4, "ENST003": 0.7})
    tx2gene = {"ENST001": "ENSG001", "ENST002": "ENSG001", "ENST003": "ENSG002"}
    gene_scores, rep_tx = project_tx_to_gene(tx_scores, tx2gene, method="max")
    assert gene_scores["ENSG001"] == pytest.approx(0.9)
    assert rep_tx["ENSG001"] == "ENST001"


def test_build_target_list_gene(gene2symbol):
    edges = pd.DataFrame(
        {
            "target_gene": ["ENSG001", "ENSG001", "ENSG002", "ENSG003"],
            "mean_importance": [0.9, 0.7, 0.8, 0.5],
        }
    )
    top_df, symbols, meta = build_target_list(
        edges, "target_gene", "mean_importance", target_type="gene",
        gene2symbol=gene2symbol, K=10
    )
    assert len(top_df) > 0
    assert len(symbols) > 0
    assert meta["n_edges"] == 4


def test_build_target_list_transcript(tx2gene, gene2symbol):
    edges = pd.DataFrame(
        {
            "target_transcript": ["ENST001", "ENST003", "ENST004"],
            "mean_importance": [0.9, 0.8, 0.7],
        }
    )
    top_df, symbols, meta = build_target_list(
        edges, "target_transcript", "mean_importance",
        target_type="transcript", tx2gene=tx2gene,
        gene2symbol=gene2symbol, K=10
    )
    assert len(symbols) > 0


def test_build_target_list_empty():
    top_df, symbols, meta = build_target_list(
        pd.DataFrame(), "target", "importance"
    )
    assert len(top_df) == 0
    assert symbols == []
