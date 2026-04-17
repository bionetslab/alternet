"""Keywords analysis utilities for g:Profiler enrichment results."""

import re

import pandas as pd


def load_gprofiler(filepath: str, dataset_id: str) -> pd.DataFrame:
    """Load a g:Profiler CSV file and return a long-form DataFrame.

    Parses the ``highlighted`` column (comma-separated per query) and melts
    the result to one row per (term, query).

    Args:
        filepath: Path to the g:Profiler CSV output file.
        dataset_id: Identifier for this dataset (added as a column).

    Returns:
        Long-form DataFrame with columns: ``dataset_id``, ``query``,
        ``source``, ``term_id``, ``term_name``, ``p_value``, ``highlighted``.
    """
    raw = pd.read_csv(filepath)

    raw["hl_L1"] = raw["highlighted"].str.split(",").str[0].str.strip() == "true"
    raw["hl_L2"] = raw["highlighted"].str.split(",").str[1].str.strip() == "true"
    raw["hl_L3"] = raw["highlighted"].str.split(",").str[2].str.strip() == "true"

    col_map = {
        "adjusted_p_value__L1_Canonical": "pval_L1",
        "adjusted_p_value__L2_Source_Isoform": "pval_L2",
        "adjusted_p_value__L3_Target_Isoform": "pval_L3",
    }
    raw = raw.rename(columns=col_map)

    rows = []
    for _, r in raw.iterrows():
        for L in ["L1", "L2", "L3"]:
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "query": L,
                    "term_id": r.get("term_id", r.get("native", "")),
                    "term_name": r["term_name"],
                    "source": r.get("source", ""),
                    "p_value": r[f"pval_{L}"],
                    "highlighted": r[f"hl_{L}"],
                }
            )
    return pd.DataFrame(rows)


def normalize_term(s: str) -> str:
    """Normalize a term name for keyword matching.

    Lowercases the string, replaces hyphens, underscores and slashes with
    spaces, and removes extra whitespace.

    Args:
        s: Term name string.

    Returns:
        Normalized string.
    """
    if pd.isna(s):
        return ""
    s = str(s).lower()
    s = re.sub(r"[-_/]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def compile_pat(keyword_list: list) -> re.Pattern:
    """Compile a regex pattern from a list of keywords.

    Short words (length <= 3) are wrapped with word boundaries.

    Args:
        keyword_list: List of keyword strings.

    Returns:
        Compiled case-insensitive regex pattern.
    """
    parts = []
    for kw in keyword_list:
        escaped = re.escape(kw.lower())
        if len(kw) <= 3:
            parts.append(r"\b" + escaped + r"\b")
        else:
            parts.append(escaped)
    return re.compile("|".join(parts), re.IGNORECASE)
