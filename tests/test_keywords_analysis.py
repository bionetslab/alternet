"""Tests for alternet.keywords_analysis module."""

import re

import pandas as pd
import pytest

from alternet.keywords_analysis import compile_pat, load_gprofiler, normalize_term


def test_normalize_term():
    assert normalize_term("Oxidative-Phosphorylation") == "oxidative phosphorylation"
    assert normalize_term("RNA_processing") == "rna processing"
    assert normalize_term("a/b") == "a b"
    assert normalize_term(None) == ""


def test_compile_pat_short_word():
    pat = compile_pat(["rna", "mrna"])
    assert pat.search("rna processing") is not None
    assert pat.search("mrna splicing") is not None
    # Short words use word boundaries - "crna" should not match "rna" boundary pattern
    assert pat.search("crna") is None


def test_compile_pat_long_word():
    pat = compile_pat(["spliceosome", "ribosome"])
    assert pat.search("ribosome biogenesis") is not None
    assert pat.search("spliceosome assembly") is not None


def test_load_gprofiler(tmp_path):
    csv_content = (
        "term_id,term_name,source,highlighted,"
        "adjusted_p_value__L1_Canonical,"
        "adjusted_p_value__L2_Source_Isoform,"
        "adjusted_p_value__L3_Target_Isoform\n"
        "GO:0001,test term,GO:BP,\"true, false, true\",0.01,0.05,0.02\n"
    )
    csv_file = tmp_path / "test_gprofiler.csv"
    csv_file.write_text(csv_content)
    df = load_gprofiler(str(csv_file), "TEST_DS")
    assert len(df) == 3  # 3 queries: L1, L2, L3
    assert set(df["query"]) == {"L1", "L2", "L3"}
    assert (df["dataset_id"] == "TEST_DS").all()
    l1_row = df[df["query"] == "L1"].iloc[0]
    assert l1_row["highlighted"] == True
    l2_row = df[df["query"] == "L2"].iloc[0]
    assert l2_row["highlighted"] == False
