"""Statistical analysis utilities for AlterNet 2.0 results."""

import os
import os.path as op

import numpy as np
import pandas as pd


def load_dataset_tables(
    dtype: str,
    group: str,
    part3_results: dict,
) -> dict:
    """Load plausibility filtering result tables for a dataset.

    Args:
        dtype: Dataset type, e.g. ``'MAGNet'`` or ``'GTEx'``.
        group: Dataset group name, e.g. ``'NF'`` or ``'Spleen'``.
        part3_results: Dict mapping dtype to results directory path.

    Returns:
        Dictionary with keys ``t1`` (set_a), ``t2`` (set_b), ``t3`` (set_c),
        each containing the corresponding DataFrame.
    """
    rdir = part3_results[dtype]
    prefix = group
    d = {}
    for tname, fname in [("t1", "set_a"), ("t2", "set_b"), ("t3", "set_c")]:
        path = op.join(rdir, f"{prefix}_{fname}_plausible.tsv")
        d[tname] = pd.read_csv(path, sep="\t")
    return d


def compute_section_A(
    data: dict,
    tx_appris_label: dict = None,
    tx_trifid: dict = None,
    gene2txs: dict = None,
    seed: int = 42,
) -> dict:
    """Compute target isoform resolution statistics from Set B.

    Args:
        data: Dict with key ``t2`` containing the Set B DataFrame.
        tx_appris_label: Optional dict mapping transcript_id to APPRIS label.
        tx_trifid: Optional dict mapping transcript_id to TRIFID score.
        gene2txs: Optional dict mapping gene_id to set of transcript IDs.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with ``stats`` and ``plot_data`` keys.
    """
    t2 = data["t2"].copy()
    st, pl = {}, {}

    st["n_edges"] = len(t2)
    n_iso = (t2["target_category"] == "target_isoform_specific").sum()
    st["n_isospec"] = int(n_iso)
    st["frac_isospec"] = n_iso / max(len(t2), 1) * 100

    iso = t2[t2["target_category"] == "target_isoform_specific"]

    if "tgt_n_isoforms" in iso.columns:
        pl["ntgt_isospec"] = (
            iso.groupby("target_gene")["tgt_n_isoforms"].first().dropna().values
        )
        st["ntgt_median_isospec"] = (
            float(np.median(pl["ntgt_isospec"])) if len(pl["ntgt_isospec"]) else np.nan
        )
    else:
        pl["ntgt_isospec"] = np.array([])
        st["ntgt_median_isospec"] = np.nan

    if "tgt_dominance" in iso.columns:
        pl["dom_isospec"] = iso["tgt_dominance"].dropna().values
        st["dom_median_isospec"] = (
            float(np.median(pl["dom_isospec"])) if len(pl["dom_isospec"]) else np.nan
        )
    else:
        pl["dom_isospec"] = np.array([])
        st["dom_median_isospec"] = np.nan

    if "target_tx_set" in t2.columns:
        iso_with_tx = iso[iso["target_tx_set"].astype(str).str.len() > 0].copy()
        if len(iso_with_tx) > 0:
            iso_with_tx["target_tx_list"] = iso_with_tx["target_tx_set"].str.split(",")
            exploded = iso_with_tx.explode("target_tx_list")
            exploded = exploded[exploded["target_tx_list"].str.strip().str.len() > 0]
            if "regulator_gene" in exploded.columns:
                regs_per_tx = exploded.groupby("target_tx_list")["regulator_gene"].nunique()
                pl["regs_per_tgt_tx"] = regs_per_tx.values
                st["regs_per_tgt_tx_median"] = float(regs_per_tx.median())

    return {"stats": st, "plot_data": pl}


def compute_section_B(data: dict) -> dict:
    """Compute SF splicing support statistics from Set C.

    Args:
        data: Dict with key ``t3`` containing the Set C DataFrame.

    Returns:
        Dictionary with ``stats`` and ``plot_data`` keys.
    """
    t3 = data["t3"].copy()
    st, pl = {}, {}

    st["n_edges"] = len(t3)
    n_sup = (t3["sf_category"] == "sf_splicing_supported_specific").sum()
    st["n_supported"] = int(n_sup)
    st["frac_supported"] = n_sup / max(len(t3), 1) * 100

    sup = t3[t3["sf_category"] == "sf_splicing_supported_specific"]

    if "tgt_n_isoforms" in sup.columns:
        pl["ntgt_sup"] = (
            sup.groupby("target_gene")["tgt_n_isoforms"].first().dropna().values
        )
        st["ntgt_median_sup"] = (
            float(np.median(pl["ntgt_sup"])) if len(pl["ntgt_sup"]) else np.nan
        )
    else:
        pl["ntgt_sup"] = np.array([])
        st["ntgt_median_sup"] = np.nan

    if "dominance" in sup.columns:
        pl["dom_sup"] = sup["dominance"].dropna().values
        st["dom_median_sup"] = (
            float(np.median(pl["dom_sup"])) if len(pl["dom_sup"]) else np.nan
        )
    else:
        pl["dom_sup"] = np.array([])
        st["dom_median_sup"] = np.nan

    if "target_tx_set" in t3.columns:
        sup_with_tx = sup[sup["target_tx_set"].astype(str).str.len() > 0].copy()
        if len(sup_with_tx) > 0:
            sup_with_tx["target_tx_list"] = sup_with_tx["target_tx_set"].str.split(",")
            exploded = sup_with_tx.explode("target_tx_list")
            exploded = exploded[exploded["target_tx_list"].str.strip().str.len() > 0]
            if "sf_gene" in exploded.columns:
                sfs_per_tx = exploded.groupby("target_tx_list")["sf_gene"].nunique()
                pl["sfs_per_tgt_tx"] = sfs_per_tx.values
                st["sfs_per_tgt_tx_median"] = float(sfs_per_tx.median())

    return {"stats": st, "plot_data": pl}


def compute_section_C(
    data: dict,
    tx_appris_label: dict,
    tx_trifid: dict,
    gene2txs: dict,
    seed: int = 42,
) -> dict:
    """Compute source annotation statistics from Set A.

    Args:
        data: Dict with key ``t1`` containing the Set A DataFrame.
        tx_appris_label: Dict mapping transcript_id to APPRIS annotation label.
        tx_trifid: Dict mapping transcript_id to TRIFID score.
        gene2txs: Dict mapping gene_id to set of transcript IDs.
        seed: Random seed for permutation tests.

    Returns:
        Dictionary with ``stats`` and ``plot_data`` keys.
    """
    rng = np.random.default_rng(seed)
    t1 = data["t1"].copy()
    st, pl = {}, {}

    def _classify_appris(label):
        if pd.isna(label) or label == "":
            return "UNKNOWN"
        s = str(label).upper()
        if "PRINCIPAL" in s:
            return "PRINCIPAL"
        if "ALTERNATIVE" in s:
            return "ALTERNATIVE"
        if "MINOR" in s:
            return "MINOR"
        return "UNKNOWN"

    t1["appris_coarse"] = t1["reg_appris"].apply(_classify_appris)
    t1["trifid_score"] = t1["best_tx"].map(tx_trifid) if tx_trifid else np.nan

    iso_t1 = t1[t1["source_category"] == "source_isoform_specific"]
    iso_tx = t1[["best_tx", "appris_coarse", "trifid_score"]].drop_duplicates(
        subset=["best_tx"]
    )
    bg_tx = t1[["best_tx", "appris_coarse", "trifid_score"]].drop_duplicates(
        subset=["best_tx"]
    )

    st["n_iso_tx"] = len(iso_tx)
    st["n_bg_tx"] = len(bg_tx)
    st["n_edges"] = len(t1)

    appris_cats = ["PRINCIPAL", "ALTERNATIVE", "MINOR"]
    for prefix_label, tx_df in [("iso", iso_tx), ("bg", bg_tx)]:
        known = tx_df[tx_df["appris_coarse"].isin(appris_cats)]
        total_known = max(len(known), 1)
        for cat in appris_cats:
            st[f"{prefix_label}_appris_{cat}_pct"] = (
                (known["appris_coarse"] == cat).sum() / total_known * 100
            )

    trifid_bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
    trifid_labels = [
        "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
        "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0",
    ]
    for prefix_label, tx_df in [("iso", iso_tx), ("bg", bg_tx)]:
        scores = tx_df["trifid_score"].dropna().values
        counts, _ = np.histogram(scores, bins=trifid_bins)
        total = max(counts.sum(), 1)
        for j, lbl in enumerate(trifid_labels):
            st[f"{prefix_label}_trifid_bin_{lbl}"] = counts[j] / total * 100

    return {"stats": st, "plot_data": pl}
