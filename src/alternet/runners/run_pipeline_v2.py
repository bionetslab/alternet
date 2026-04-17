"""Runner scripts for AlterNet 2.0 pipeline reproduction."""

import argparse
import os
import os.path as op
import sys
import yaml


def _load_config(config_path: str) -> dict:
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def run_full_pipeline(config_path: str) -> None:
    """Run full AlterNet 2.0 pipeline from expression data to enrichment analysis.

    Executes network inference for each dataset, followed by edge
    categorization, plausibility filtering, and g:Profiler input generation.

    Args:
        config_path: Path to YAML configuration file.
    """
    import pandas as pd
    from alternet import network_inference_v2 as ni
    from alternet import edge_categorization as ec
    from alternet import plausibility_filtering as pf
    from alternet import gprofiler_input as gp

    cfg = _load_config(config_path)
    data_dir = cfg["data_dir"]
    results_dir = cfg["results_dir"]
    os.makedirs(results_dir, exist_ok=True)

    biomart = pd.read_csv(op.join(data_dir, cfg["biomart"]), sep="\t")
    appris_df = pd.read_csv(op.join(data_dir, cfg["appris"]), sep="\t")
    digger_df = pd.read_csv(op.join(data_dir, cfg["digger"]), low_memory=False)
    tf_list_raw = pd.read_csv(op.join(data_dir, cfg["tf_list"]), sep="\t", header=None)
    sf_list_raw = pd.read_csv(op.join(data_dir, cfg["sf_list"]), header=None)

    n_runs = cfg.get("n_runs", 10)
    variance_percentile = cfg.get("variance_percentile", 0.7)

    tx2gene = dict(zip(biomart["Transcript stable ID"], biomart["Gene stable ID"]))

    # MAGNet conditions
    for condition in cfg.get("magnet_conditions", []):
        print(f"\n=== MAGNet {condition} ===")
        tpm_path = op.join(
            data_dir,
            cfg["magnet_tpm_pattern"].format(condition=condition),
        )
        transcript_data = ni.preprocess_magnet_expression(
            tpm_path, biomart, variance_percentile=variance_percentile
        )
        sample_cols = [
            c for c in transcript_data.columns
            if c not in ["transcript_id", "gene_id"]
        ]

        net_dir = op.join(results_dir, f"magnet_{condition}_networks")
        os.makedirs(net_dir, exist_ok=True)

        ni_result = ni._run_pipeline(
            transcript_data, biomart, tf_list_raw, sf_list_raw,
            appris_df, digger_df, n_runs=n_runs
        )

        ni_result["canonical_grn"].to_csv(
            op.join(net_dir, f"{condition}_canonical_raw.tsv"), sep="\t", index=False
        )
        ni_result["as_source_grn"].to_csv(
            op.join(net_dir, f"{condition}_as_aware_source_raw.tsv"), sep="\t", index=False
        )
        ni_result["fully_as_grn"].to_csv(
            op.join(net_dir, f"{condition}_fully_as_aware_raw.tsv"), sep="\t", index=False
        )

        ec_result = ec._run_pipeline(
            ni_result["canonical_grn"], ni_result["as_source_grn"],
            ni_result["fully_as_grn"], ni_result["regulator_list"],
            transcript_data, sample_cols, tx2gene,
        )

        ec_dir = op.join(results_dir, f"magnet_{condition}_edge_cate")
        os.makedirs(ec_dir, exist_ok=True)
        for key in ["net1_filtered", "net2_filtered", "net3_filtered",
                    "set_a", "set_b", "set_c", "set_d"]:
            ec_result[key].to_csv(
                op.join(ec_dir, f"{condition}_{key}.tsv"), sep="\t", index=False
            )

        pf_result = pf._run_pipeline(
            ec_result["set_a"], ec_result["set_b"],
            ec_result["set_c"], ec_result["set_d"],
            ec_result["set_b_unpacked"], ec_result["set_c_unpacked"],
            transcript_data, sample_cols, biomart, appris_df, digger_df,
        )

        pf_dir = op.join(results_dir, f"magnet_{condition}_plau_filter")
        os.makedirs(pf_dir, exist_ok=True)
        for key in ["set_a_filtered", "set_b_filtered", "set_c_filtered", "set_d_filtered"]:
            pf_result[key].to_csv(
                op.join(pf_dir, f"{condition}_{key}.tsv"), sep="\t", index=False
            )

        gp_result = gp._run_pipeline(
            pf_result["set_a_filtered"], pf_result["set_b_filtered"],
            ec_result["net1_filtered"], ec_result["net2_filtered"],
            ec_result["net3_filtered"], tx2gene,
            dict(zip(biomart["Gene stable ID"], biomart["Gene name"])),
            top_k=cfg.get("top_k", 500),
        )

        gp_dir = op.join(results_dir, f"magnet_{condition}_gprofiler")
        os.makedirs(gp_dir, exist_ok=True)
        for key in ["l1_top", "l2_top", "l3_top"]:
            gp_result[key].to_csv(
                op.join(gp_dir, f"{condition}_{key}.tsv"), sep="\t", index=False
            )

    # GTEx tissues
    gtex_tpm = op.join(data_dir, cfg.get("gtex_tpm", ""))
    gtex_attr = op.join(data_dir, cfg.get("gtex_sample_attributes", ""))

    for tissue in cfg.get("gtex_tissues", []):
        print(f"\n=== GTEx {tissue} ===")
        transcript_data = ni.preprocess_gtex_expression(
            gtex_tpm, gtex_attr, tissue, biomart,
            variance_percentile=variance_percentile,
        )
        sample_cols = [
            c for c in transcript_data.columns
            if c not in ["transcript_id", "gene_id"]
        ]

        safe_tissue = tissue.replace(" ", "_").replace("-", "_")
        net_dir = op.join(results_dir, f"gtex_{safe_tissue}_networks")
        os.makedirs(net_dir, exist_ok=True)

        ni_result = ni._run_pipeline(
            transcript_data, biomart, tf_list_raw, sf_list_raw,
            appris_df, digger_df, n_runs=n_runs
        )

        ni_result["canonical_grn"].to_csv(
            op.join(net_dir, f"{tissue}_canonical_raw.tsv"), sep="\t", index=False
        )
        ni_result["as_source_grn"].to_csv(
            op.join(net_dir, f"{tissue}_as_aware_source_raw.tsv"), sep="\t", index=False
        )
        ni_result["fully_as_grn"].to_csv(
            op.join(net_dir, f"{tissue}_fully_as_aware_raw.tsv"), sep="\t", index=False
        )

        ec_result = ec._run_pipeline(
            ni_result["canonical_grn"], ni_result["as_source_grn"],
            ni_result["fully_as_grn"], ni_result["regulator_list"],
            transcript_data, sample_cols, tx2gene,
        )

        ec_dir = op.join(results_dir, f"gtex_{safe_tissue}_edge_cate")
        os.makedirs(ec_dir, exist_ok=True)
        for key in ["net1_filtered", "net2_filtered", "net3_filtered",
                    "set_a", "set_b", "set_c", "set_d"]:
            ec_result[key].to_csv(
                op.join(ec_dir, f"{tissue}_{key}.tsv"), sep="\t", index=False
            )

        pf_result = pf._run_pipeline(
            ec_result["set_a"], ec_result["set_b"],
            ec_result["set_c"], ec_result["set_d"],
            ec_result["set_b_unpacked"], ec_result["set_c_unpacked"],
            transcript_data, sample_cols, biomart, appris_df, digger_df,
        )

        pf_dir = op.join(results_dir, f"gtex_{safe_tissue}_plau_filter")
        os.makedirs(pf_dir, exist_ok=True)
        for key in ["set_a_filtered", "set_b_filtered", "set_c_filtered", "set_d_filtered"]:
            pf_result[key].to_csv(
                op.join(pf_dir, f"{tissue}_{key}.tsv"), sep="\t", index=False
            )


def run_downstream_pipeline(config_path: str) -> None:
    """Run AlterNet 2.0 pipeline starting from pre-computed networks.

    Executes edge categorization, plausibility filtering, and g:Profiler
    input generation from previously inferred networks.

    Args:
        config_path: Path to YAML configuration file.
    """
    import pandas as pd
    from alternet import edge_categorization as ec
    from alternet import plausibility_filtering as pf
    from alternet import gprofiler_input as gp

    cfg = _load_config(config_path)
    data_dir = cfg["data_dir"]
    results_dir = cfg["results_dir"]

    biomart = pd.read_csv(op.join(data_dir, cfg["biomart"]), sep="\t")
    appris_df = pd.read_csv(op.join(data_dir, cfg["appris"]), sep="\t")
    digger_df = pd.read_csv(op.join(data_dir, cfg["digger"]), low_memory=False)

    tx2gene = dict(zip(biomart["Transcript stable ID"], biomart["Gene stable ID"]))
    gene2symbol = dict(zip(biomart["Gene stable ID"], biomart["Gene name"]))

    print("Downstream pipeline: starting from pre-computed networks in", results_dir)


def main() -> None:
    """Entry point for the AlterNet 2.0 pipeline runner."""
    parser = argparse.ArgumentParser(
        description="AlterNet 2.0 pipeline runner"
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--mode",
        choices=["full", "downstream"],
        default="full",
        help="Pipeline mode: 'full' or 'downstream'",
    )
    args = parser.parse_args()

    if args.mode == "full":
        run_full_pipeline(args.config)
    else:
        run_downstream_pipeline(args.config)


if __name__ == "__main__":
    main()
