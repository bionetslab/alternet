"""Tests for alternet.splicing_factors module."""

import json
import os

import pandas as pd
import pytest

from alternet.splicing_factors import (
    build_sf_database,
    load_encode_rbp,
    load_kegg_spliceosome,
    load_spliceaid_f,
    load_splicinglore,
)


# ---- helpers ----------------------------------------------------------------

def _write_tmp(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(content)


# ---- individual loaders -----------------------------------------------------

def test_load_spliceaid_f(tmp_path):
    tsv = tmp_path / "saf.tsv"
    tsv.write_text("Gene\tOther\nSF1\tval1\nSF2\tval2\n")
    result = load_spliceaid_f(str(tsv))
    assert result == ["SF1", "SF2"]


def test_load_kegg_spliceosome(tmp_path):
    data = {"KEGG_SPLICEOSOME": {"geneSymbols": ["GENE_A", "GENE_B"]}}
    jf = tmp_path / "kegg.json"
    jf.write_text(json.dumps(data))
    result = load_kegg_spliceosome(str(jf))
    assert result == ["GENE_A", "GENE_B"]


# ---- build_sf_database (mocked loaders) ----------------------------------------

def test_build_sf_database(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "alternet.splicing_factors.load_spliceaid_f", lambda p: ["SF1", "SF2"]
    )
    monkeypatch.setattr(
        "alternet.splicing_factors.load_splicinglore", lambda p, **kw: ["SF2", "SF3"]
    )
    monkeypatch.setattr(
        "alternet.splicing_factors.load_kegg_spliceosome", lambda p: ["SF3", "SF4"]
    )
    monkeypatch.setattr(
        "alternet.splicing_factors.load_encode_rbp", lambda p: ["SF1", "SF4", "SF5"]
    )

    df = build_sf_database("a", "b", "c", "d")

    assert set(df["Splicing_Factor"]) == {"SF1", "SF2", "SF3", "SF4", "SF5"}
    assert "Source_Count" in df.columns
    assert df.loc[df["Splicing_Factor"] == "SF4", "Source_Count"].values[0] == 2
