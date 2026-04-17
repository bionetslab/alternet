"""Network inference pipeline for AlterNet 2.0 (TF + SF regulators)."""

import time
import numpy as np
import pandas as pd


def map_sf_ids(sf_list_raw: pd.DataFrame, biomart: pd.DataFrame) -> pd.DataFrame:
    """Map SF gene names to Ensembl gene and transcript IDs.

    Args:
        sf_list_raw: DataFrame with a single column of SF gene names.
        biomart: BioMart reference DataFrame.

    Returns:
        DataFrame with columns ``SF``, ``Gene stable ID``, ``Transcript stable ID``.
    """
    sf_list_raw = sf_list_raw.copy()
    sf_list_raw.columns = ["SF"]
    sf_list = sf_list_raw.merge(biomart, left_on="SF", right_on="Gene name")
    sf_list = sf_list.loc[
        :, ["SF", "Gene stable ID", "Transcript stable ID"]
    ].drop_duplicates()
    return sf_list


def combine_tf_sf_lists(
    tf_list: pd.DataFrame, sf_list: pd.DataFrame
) -> pd.DataFrame:
    """Combine TF and SF lists, marking overlapping genes as TF_SF.

    Args:
        tf_list: TF DataFrame with columns ``TF``, ``Gene stable ID``,
            ``Transcript stable ID``.
        sf_list: SF DataFrame with columns ``SF``, ``Gene stable ID``,
            ``Transcript stable ID``.

    Returns:
        Combined DataFrame with columns ``Regulator_name``, ``Gene stable ID``,
        ``Transcript stable ID``, ``Regulator_type``.
    """
    tf_list = tf_list.copy()
    sf_list = sf_list.copy()

    tf_list["Regulator_type"] = "TF"
    sf_list["Regulator_type"] = "SF"
    tf_list = tf_list.rename(columns={"TF": "Regulator_name"})
    sf_list = sf_list.rename(columns={"SF": "Regulator_name"})

    tf_genes = set(tf_list["Gene stable ID"])
    sf_genes = set(sf_list["Gene stable ID"])
    overlap_genes = tf_genes & sf_genes

    combined = pd.concat([tf_list, sf_list], ignore_index=True)
    combined.loc[
        combined["Gene stable ID"].isin(overlap_genes), "Regulator_type"
    ] = "TF_SF"
    combined = combined.drop_duplicates(subset=["Transcript stable ID"], keep="first")
    return combined


def preprocess_gtex_expression(
    gtex_tpm_path: str,
    sample_attr_path: str,
    tissue: str,
    biomart: pd.DataFrame,
    variance_percentile: float = 0.7,
) -> pd.DataFrame:
    """Load and filter GTEx expression data for a given tissue.

    Filters to tissue samples, protein-coding transcripts, removes low-expression
    transcripts, and applies variance filtering.

    Args:
        gtex_tpm_path: Path to GTEx transcript TPM file.
        sample_attr_path: Path to GTEx sample attributes file.
        tissue: Tissue name as in the ``SMTSD`` column.
        biomart: BioMart reference DataFrame.
        variance_percentile: Quantile below which transcripts are filtered out.

    Returns:
        DataFrame with columns ``transcript_id``, ``gene_id``, and sample columns.
    """
    sample_attributes = pd.read_csv(sample_attr_path, sep="\t", low_memory=False)
    tissue_samples = sample_attributes[
        sample_attributes["SMTSD"] == tissue
    ]["SAMPID"].tolist()

    header_df = pd.read_csv(gtex_tpm_path, sep="\t", nrows=0)
    all_columns = header_df.columns.tolist()
    id_col_1, id_col_2 = all_columns[0], all_columns[1]
    sample_columns_in_file = all_columns[2:]
    matching_samples = [c for c in sample_columns_in_file if c in tissue_samples]

    columns_to_read = [id_col_1, id_col_2] + matching_samples
    transcript_data = pd.read_csv(
        gtex_tpm_path, sep="\t", usecols=columns_to_read
    )
    transcript_data = transcript_data.rename(
        columns={id_col_1: "transcript_id", id_col_2: "gene_id"}
    )
    transcript_data["transcript_id"] = (
        transcript_data["transcript_id"].str.split(".").str[0]
    )
    transcript_data["gene_id"] = (
        transcript_data["gene_id"].str.split(".").str[0]
    )

    sample_cols = [
        c for c in transcript_data.columns if c not in ["transcript_id", "gene_id"]
    ]

    protein_coding = biomart[
        biomart["Gene type"] == "protein_coding"
    ]["Transcript stable ID"].unique()
    transcript_data = transcript_data[
        transcript_data["transcript_id"].isin(protein_coding)
    ].copy()

    threshold = len(sample_cols) * 0.1
    zero_counts = (transcript_data[sample_cols] == 0).sum(axis=1)
    transcript_data = transcript_data[zero_counts < threshold].copy()

    expression_values = transcript_data[sample_cols].values
    log_expr = np.log1p(expression_values)
    variances = np.var(log_expr, axis=1)
    variance_threshold = np.quantile(variances, variance_percentile)
    transcript_data = transcript_data[variances > variance_threshold].copy()

    return transcript_data


def preprocess_magnet_expression(
    magnet_tpm_path: str,
    biomart: pd.DataFrame,
    apply_variance_filter: bool = True,
    variance_percentile: float = 0.7,
) -> pd.DataFrame:
    """Load and filter MAGNet expression data.

    Applies protein-coding filtering and optional variance filtering.

    Args:
        magnet_tpm_path: Path to MAGNet pre-filtered TPM TSV file.
        biomart: BioMart reference DataFrame.
        apply_variance_filter: Whether to apply variance filtering.
        variance_percentile: Quantile below which transcripts are filtered out.

    Returns:
        DataFrame with columns ``transcript_id``, ``gene_id``, and sample columns.
    """
    transcript_data = pd.read_csv(magnet_tpm_path, sep="\t")

    sample_cols = [
        c for c in transcript_data.columns if c not in ["transcript_id", "gene_id"]
    ]

    protein_coding = biomart[
        biomart["Gene type"] == "protein_coding"
    ]["Transcript stable ID"].unique()
    transcript_data = transcript_data[
        transcript_data["transcript_id"].isin(protein_coding)
    ].copy()

    if apply_variance_filter:
        expression_values = transcript_data[sample_cols].values
        log_expr = np.log1p(expression_values)
        variances = np.var(log_expr, axis=1)
        variance_threshold = np.quantile(variances, variance_percentile)
        transcript_data = transcript_data[variances > variance_threshold].copy()

    return transcript_data


def _run_pipeline(
    transcript_data: pd.DataFrame,
    biomart: pd.DataFrame,
    tf_list_raw: pd.DataFrame,
    sf_list_raw: pd.DataFrame,
    appris_df: pd.DataFrame,
    digger_df: pd.DataFrame,
    n_runs: int = 10,
) -> dict:
    """Run network inference for all three networks.

    Args:
        transcript_data: Preprocessed transcript TPM data with ``transcript_id``
            and ``gene_id`` columns plus sample columns.
        biomart: BioMart reference DataFrame.
        tf_list_raw: Raw TF list DataFrame (single column of gene names).
        sf_list_raw: Raw SF list DataFrame (single column of gene names).
        appris_df: APPRIS annotation DataFrame.
        digger_df: DIGGER annotation DataFrame.
        n_runs: Number of GRNBoost2 runs.

    Returns:
        Dictionary with keys: ``canonical_grn``, ``as_source_grn``,
        ``fully_as_grn``, ``tf_list``, ``sf_list``, ``regulator_list``,
        ``tf_isoform_categories``, ``regulator_isoform_categories``,
        ``target_isoform_categories``, ``tf_gene_categories``,
        ``regulator_gene_categories``, ``target_gene_categories``,
        ``summary_stats``.
    """
    from alternet.annotation import map_tf_ids
    import alternet.postprocessing as postprocessing
    from alternet.data_preprocessing import create_hybrid_data, standardize_dataframe
    from alternet.inference import inference

    sample_cols = [
        c for c in transcript_data.columns if c not in ["transcript_id", "gene_id"]
    ]

    tf_list = map_tf_ids(tf_list_raw, biomart)
    sf_list = map_sf_ids(sf_list_raw, biomart)
    regulator_list = combine_tf_sf_lists(tf_list, sf_list)

    gene_data = (
        transcript_data.groupby("gene_id")[sample_cols].sum().reset_index()
    )

    gene_data_matrix = gene_data.set_index("gene_id")[sample_cols].T
    transcript_data_matrix = transcript_data.set_index("transcript_id")[sample_cols].T

    gene_data_scaled = standardize_dataframe(gene_data_matrix)
    transcript_data_scaled = standardize_dataframe(transcript_data_matrix)

    nan_cols = transcript_data_scaled.columns[
        transcript_data_scaled.isna().any()
    ].tolist()
    zero_var_cols = transcript_data_matrix.columns[
        transcript_data_matrix.std() == 0
    ].tolist()
    bad_transcripts = set(nan_cols + zero_var_cols)
    if bad_transcripts:
        good_transcripts = [
            c for c in transcript_data_scaled.columns if c not in bad_transcripts
        ]
        transcript_data_scaled = transcript_data_scaled[good_transcripts]
        transcript_data_matrix = transcript_data_matrix[good_transcripts]

    transcript_data_scaled = transcript_data_scaled.fillna(0)
    gene_data_scaled = gene_data_scaled.fillna(0)

    tf_genes_in_data = list(
        set(tf_list["Gene stable ID"]) & set(gene_data_scaled.columns)
    )
    tf_transcripts_in_data = list(
        set(tf_list["Transcript stable ID"]) & set(transcript_data_scaled.columns)
    )
    regulator_transcripts_in_data = list(
        set(regulator_list["Transcript stable ID"])
        & set(transcript_data_scaled.columns)
    )
    target_genes = list(gene_data_scaled.columns)
    target_transcripts = list(transcript_data_scaled.columns)

    tf_only_transcripts = (
        set(regulator_list[regulator_list["Regulator_type"] == "TF"]["Transcript stable ID"])
        & set(transcript_data_scaled.columns)
    )
    sf_only_transcripts = (
        set(regulator_list[regulator_list["Regulator_type"] == "SF"]["Transcript stable ID"])
        & set(transcript_data_scaled.columns)
    )
    tfsf_transcripts = (
        set(regulator_list[regulator_list["Regulator_type"] == "TF_SF"]["Transcript stable ID"])
        & set(transcript_data_scaled.columns)
    )

    tf_isoform_categories = postprocessing.isoform_categorization(
        transcript_data_matrix, gene_data_matrix, tf_list
    )
    tf_gene_categories = postprocessing.get_gene_cases(tf_isoform_categories)

    regulator_isoform_categories = postprocessing.isoform_categorization(
        transcript_data_matrix, gene_data_matrix, regulator_list
    )
    regulator_gene_categories = postprocessing.get_gene_cases(regulator_isoform_categories)

    target_list = biomart[
        biomart["Transcript stable ID"].isin(transcript_data_scaled.columns)
    ][["Gene stable ID", "Transcript stable ID"]].drop_duplicates()
    target_list = target_list[
        target_list["Gene stable ID"].isin(gene_data_scaled.columns)
    ]
    target_isoform_categories = postprocessing.isoform_categorization(
        transcript_data_matrix, gene_data_matrix, target_list
    )
    target_gene_categories = postprocessing.get_gene_cases(target_isoform_categories)

    runtime = {}

    start = time.monotonic()
    canonical_grn = inference(
        gene_data=gene_data_scaled,
        tf_list=tf_genes_in_data,
        target_names="all",
        n_runs=n_runs,
    )
    runtime["canonical"] = time.monotonic() - start

    hybrid_data = create_hybrid_data(
        transcript_data_matrix, gene_data_matrix, tf_list
    )
    start = time.monotonic()
    as_source_grn = inference(
        gene_data=hybrid_data,
        tf_list=tf_transcripts_in_data,
        target_names=target_genes,
        n_runs=n_runs,
    )
    runtime["as_aware_source"] = time.monotonic() - start

    start = time.monotonic()
    fully_as_grn = inference(
        gene_data=transcript_data_scaled,
        tf_list=regulator_transcripts_in_data,
        target_names="all",
        n_runs=n_runs,
    )
    runtime["fully_as_aware"] = time.monotonic() - start

    runtime["total"] = (
        runtime["canonical"] + runtime["as_aware_source"] + runtime["fully_as_aware"]
    )

    summary_stats = {
        "n_samples": len(sample_cols),
        "n_genes": len(gene_data_scaled.columns),
        "n_transcripts": len(transcript_data_scaled.columns),
        "n_tf_genes": len(tf_genes_in_data),
        "n_tf_transcripts": len(tf_transcripts_in_data),
        "n_sf_transcripts": len(sf_only_transcripts),
        "n_tfsf_transcripts": len(tfsf_transcripts),
        "n_regulator_transcripts": len(regulator_transcripts_in_data),
        "network1_edges": len(canonical_grn),
        "network2_edges": len(as_source_grn),
        "network3_edges": len(fully_as_grn),
        "runtime_canonical_min": runtime["canonical"] / 60,
        "runtime_as_source_min": runtime["as_aware_source"] / 60,
        "runtime_fully_as_min": runtime["fully_as_aware"] / 60,
        "runtime_total_min": runtime["total"] / 60,
    }

    return {
        "canonical_grn": canonical_grn,
        "as_source_grn": as_source_grn,
        "fully_as_grn": fully_as_grn,
        "tf_list": tf_list,
        "sf_list": sf_list,
        "regulator_list": regulator_list,
        "tf_isoform_categories": tf_isoform_categories,
        "regulator_isoform_categories": regulator_isoform_categories,
        "target_isoform_categories": target_isoform_categories,
        "tf_gene_categories": tf_gene_categories,
        "regulator_gene_categories": regulator_gene_categories,
        "target_gene_categories": target_gene_categories,
        "summary_stats": summary_stats,
    }
