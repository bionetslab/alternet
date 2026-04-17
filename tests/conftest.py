"""Shared fixtures for alternet tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def biomart_df():
    """Minimal BioMart DataFrame."""
    return pd.DataFrame(
        {
            "Transcript stable ID": [
                "ENST001", "ENST002", "ENST003", "ENST004", "ENST005",
            ],
            "Gene stable ID": [
                "ENSG001", "ENSG001", "ENSG002", "ENSG003", "ENSG003",
            ],
            "Gene name": ["GENE1", "GENE1", "GENE2", "GENE3", "GENE3"],
            "Gene type": [
                "protein_coding", "protein_coding", "protein_coding",
                "protein_coding", "protein_coding",
            ],
        }
    )


@pytest.fixture
def transcript_tpm_df():
    """Minimal transcript TPM DataFrame with 30 samples."""
    rng = np.random.default_rng(42)
    n_tx = 5
    n_samples = 30
    sample_cols = [f"SAMPLE{i:02d}" for i in range(n_samples)]
    data = {"transcript_id": ["ENST001", "ENST002", "ENST003", "ENST004", "ENST005"],
            "gene_id": ["ENSG001", "ENSG001", "ENSG002", "ENSG003", "ENSG003"]}
    for col in sample_cols:
        data[col] = rng.exponential(scale=5.0, size=n_tx)
    return pd.DataFrame(data)


@pytest.fixture
def sample_cols(transcript_tpm_df):
    """List of sample column names."""
    return [c for c in transcript_tpm_df.columns if c not in ["transcript_id", "gene_id"]]


@pytest.fixture
def tx2gene():
    """Transcript-to-gene mapping."""
    return {
        "ENST001": "ENSG001",
        "ENST002": "ENSG001",
        "ENST003": "ENSG002",
        "ENST004": "ENSG003",
        "ENST005": "ENSG003",
    }


@pytest.fixture
def gene2symbol():
    """Gene-to-symbol mapping."""
    return {
        "ENSG001": "GENE1",
        "ENSG002": "GENE2",
        "ENSG003": "GENE3",
    }


@pytest.fixture
def minimal_net1():
    """Minimal Network 1 edge DataFrame."""
    return pd.DataFrame(
        {
            "source": ["ENSG001", "ENSG001", "ENSG002"],
            "target": ["ENSG002", "ENSG003", "ENSG003"],
            "importance": [0.8, 0.6, 0.5],
            "mean_importance": [0.8, 0.6, 0.5],
            "median_importance": [0.75, 0.55, 0.45],
            "frequency": [15, 12, 11],
        }
    )


@pytest.fixture
def minimal_net2():
    """Minimal Network 2 edge DataFrame."""
    return pd.DataFrame(
        {
            "source": ["ENST001", "ENST001", "ENST003"],
            "target": ["ENSG002", "ENSG003", "ENSG003"],
            "importance": [0.9, 0.7, 0.5],
            "mean_importance": [0.9, 0.7, 0.5],
            "median_importance": [0.85, 0.65, 0.45],
            "frequency": [15, 13, 11],
        }
    )


@pytest.fixture
def minimal_net3():
    """Minimal Network 3 edge DataFrame."""
    return pd.DataFrame(
        {
            "source": ["ENST001", "ENST001"],
            "target": ["ENST003", "ENST004"],
            "importance": [1.0, 0.8],
            "mean_importance": [1.0, 0.8],
            "median_importance": [0.95, 0.75],
            "frequency": [15, 12],
        }
    )


@pytest.fixture
def regulator_list_df():
    """Minimal regulator list DataFrame."""
    return pd.DataFrame(
        {
            "Regulator_name": ["TF1", "SF1", "TFSF1"],
            "Gene stable ID": ["ENSG001", "ENSG003", "ENSG002"],
            "Transcript stable ID": ["ENST001", "ENST004", "ENST003"],
            "Regulator_type": ["TF", "SF", "TF_SF"],
        }
    )


@pytest.fixture
def appris_df():
    """Minimal APPRIS annotation DataFrame."""
    return pd.DataFrame(
        {
            "Transcript ID": ["ENST001.1", "ENST002.1", "ENST003.1", "ENST004.1", "ENST005.1"],
            "APPRIS Annotation": [
                "PRINCIPAL:1", "ALTERNATIVE:1", "PRINCIPAL:1",
                "PRINCIPAL:1", "ALTERNATIVE:2",
            ],
            "Trifid Score": [0.9, 0.4, 0.85, 0.75, 0.3],
        }
    )


@pytest.fixture
def digger_df():
    """Minimal DIGGER annotation DataFrame."""
    return pd.DataFrame(
        {
            "Transcript stable ID": ["ENST001.1", "ENST001.1", "ENST002.1", "ENST003.1"],
            "Gene stable ID": ["ENSG001.1", "ENSG001.1", "ENSG001.1", "ENSG002.1"],
            "Exon rank in transcript": [1, 2, 1, 1],
            "Pfam ID": ["PF00001", "PF00002", None, "PF00003"],
            "Exon stable ID": ["ENSE001", "ENSE002", "ENSE003", "ENSE004"],
        }
    )
