"""Tests for alternet.network_inference_v2 module."""

import numpy as np
import pandas as pd
import pytest

from alternet.network_inference_v2 import (
    combine_tf_sf_lists,
    map_sf_ids,
)


def test_map_sf_ids(biomart_df):
    sf_raw = pd.DataFrame({"SF": ["GENE1", "GENE3"]})
    result = map_sf_ids(sf_raw, biomart_df)
    assert "SF" in result.columns
    assert "Gene stable ID" in result.columns
    assert "Transcript stable ID" in result.columns
    assert set(result["SF"]) == {"GENE1", "GENE3"}


def test_combine_tf_sf_lists(biomart_df):
    tf_list = pd.DataFrame(
        {
            "TF": ["GENE1", "GENE2"],
            "Gene stable ID": ["ENSG001", "ENSG002"],
            "Transcript stable ID": ["ENST001", "ENST003"],
        }
    )
    sf_list = pd.DataFrame(
        {
            "SF": ["GENE2", "GENE3"],
            "Gene stable ID": ["ENSG002", "ENSG003"],
            "Transcript stable ID": ["ENST003", "ENST004"],
        }
    )
    result = combine_tf_sf_lists(tf_list, sf_list)
    assert "Regulator_type" in result.columns
    assert "TF_SF" in result["Regulator_type"].values
    overlap = result[result["Gene stable ID"] == "ENSG002"]
    assert (overlap["Regulator_type"] == "TF_SF").all()
