"""Splicing factor database construction from multiple sources."""

import json
import pandas as pd


def load_spliceaid_f(path: str) -> list:
    """Read SpliceAid-F TSV file and return list of gene names.

    Args:
        path: Path to SpliceAid-F TSV file with a ``Gene`` column.

    Returns:
        List of gene names.
    """
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    return df["Gene"].dropna().tolist()


def load_splicinglore(path: str, sheet_name: str = "SF") -> list:
    """Read SplicingLore Excel file and return list of splicing factor names.

    Args:
        path: Path to SplicingLore Excel file.
        sheet_name: Sheet to read (default ``'SF'``).

    Returns:
        List of splicing factor names.
    """
    df = pd.read_excel(path, sheet_name=sheet_name)
    return df["Splicing factor"].dropna().tolist()


def load_kegg_spliceosome(path: str) -> list:
    """Read KEGG spliceosome JSON file and return list of gene symbols.

    Args:
        path: Path to KEGG JSON file.

    Returns:
        List of gene symbols.
    """
    with open(path, "r") as f:
        data = json.load(f)
    return data["KEGG_SPLICEOSOME"]["geneSymbols"]


def load_encode_rbp(path: str) -> list:
    """Read ENCODE RBP Excel file and return splicing-related gene names.

    Reads the file with ``header=1``, renames the first column to ``RBP_Name``,
    and returns genes where either ``Splicing regulation`` or ``Spliceosome``
    equals 1.

    Args:
        path: Path to ENCODE Excel file.

    Returns:
        List of gene names with splicing activity.
    """
    raw = pd.read_excel(path, header=1)
    old_name = raw.columns[0]
    raw = raw.rename(columns={old_name: "RBP_Name"})
    filtered = raw[
        (raw["Splicing regulation"].fillna(0) == 1)
        | (raw["Spliceosome"].fillna(0) == 1)
    ]
    return filtered["RBP_Name"].dropna().tolist()


def build_sf_database(
    saf_path: str,
    splicinglore_path: str,
    kegg_path: str,
    encode_path: str,
) -> pd.DataFrame:
    """Combine all four splicing factor sources into a single DataFrame.

    Args:
        saf_path: Path to SpliceAid-F TSV file.
        splicinglore_path: Path to SplicingLore Excel file.
        kegg_path: Path to KEGG spliceosome JSON file.
        encode_path: Path to ENCODE RBP Excel file.

    Returns:
        DataFrame with columns: ``Splicing_Factor``, ``In_SpliceAidF``,
        ``In_SplicingLore``, ``In_KEGG``, ``In_ENCODE``, ``Source_Count``.
    """
    saf = set(load_spliceaid_f(saf_path))
    sl = set(load_splicinglore(splicinglore_path))
    ks = set(load_kegg_spliceosome(kegg_path))
    en = set(load_encode_rbp(encode_path))

    all_sfs = sorted(saf | sl | ks | en)

    df = pd.DataFrame(all_sfs, columns=["Splicing_Factor"])
    df = df.sort_values("Splicing_Factor").reset_index(drop=True)

    df["In_SpliceAidF"] = df["Splicing_Factor"].isin(saf)
    df["In_SplicingLore"] = df["Splicing_Factor"].isin(sl)
    df["In_KEGG"] = df["Splicing_Factor"].isin(ks)
    df["In_ENCODE"] = df["Splicing_Factor"].isin(en)

    annotation_cols = ["In_SpliceAidF", "In_SplicingLore", "In_KEGG", "In_ENCODE"]
    df["Source_Count"] = df[annotation_cols].sum(axis=1)

    return df
