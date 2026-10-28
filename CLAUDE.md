# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

AlterNet infers alternative-splicing-aware gene regulatory networks (GRNs) with GRNBoost2. This version (branch `combo_net`) infers **three** networks at different resolutions, then sorts their edges by which resolution explains each interaction best. It also separates transcription factor (TF) regulators from splice factor (SF) regulators. Preprint: bioRxiv 2025.11.21.689745.

## Environment and commands

- **pixi** (linux-64 only, config in `pyproject.toml` under `[tool.pixi.*]`): installs `alternet` as an editable package.
  - `pixi run install-kernel` / `pixi run uninstall-kernel` (Jupyter kernel `pixi_alternet`)
- **pip/conda**: `pip install -e .` (hatchling, Python >= 3.10).
- `signifikante` (provides `signifikante.algo.signifikante_fdr` / `grnboost2`) is installed from the `signifikante-bootstrap` branch of `git+https://github.com/bionetslab/SignifiKANTE.git`. That branch adds `bootstraps=` (list of bootstrapped expression matrices, same shape/columns as the input) and the output columns `shap_mean`, `shap_median`, `shap_std`, `max_rank_diff`.
- `pytximport` is pinned to 0.12 (0.13 needs numpy >= 2, signifikante needs numpy < 2). It reads Salmon bootstraps.
- `splicefactor_evidence.py` imports `statsmodels` and `joblib`, but neither is a declared dependency in `pyproject.toml`.
- **Stale pixi tasks**: `pytest` (there is no `tests/` directory), `reproduce-article-results*` (`alternet.runners.run_pipeline_v2` and `configs/article_results_config.yaml` don't exist) and `construct-sf-database` (`alternet.runners.create_sf_database` doesn't exist) all fail. There is no lint config.
- The README still documents the old `alternet.runners.run_alternet.alternet_pipeline` and `tutorials/alternet_minimal.ipynb`. Both moved to `src/alternet/deprecated/`.

## Running the current pipeline (v3)

There are two CLI scripts in `src/alternet/runners/`. Both take explicit file paths.

1. **Inference**: `python -m alternet.runners.infer_networks_v3 --transcript_file ... --regulator_file data/regulator_db.tsv --biomart_file data/biomart.txt --results_path ... --experiment_name ... --canonical <out.tsv> --source_as <out.tsv> --fully_as <out.tsv>`. A network is computed only if its output file doesn't exist yet, so delete a file to recompute it.
2. **Downstream**: `python -m alternet.runners.alternet_pipeline_v3` with the same three network paths plus `--appris`, `--digger`, `--regulator_file` and `--biomart_file`. It writes `{experiment_name}.tf.tsv`, `.sf.tsv` and `.ambi.tsv` to `--results_path`.

`runners/debug_alternet.ipynb` steps through the downstream script cell by cell, using hardcoded cluster paths (`/data/bionets/...`). The `*_v2.py` runners are the previous iteration, with two networks and p-value filtering. `run_mgnet.sh` points at a script that no longer exists.

### Single-network runners (signifikante P-values)

All three call `inference.compute_grn` (signifikante_fdr, medoid/kmeans, BH-adjusted `pvalue_bh`) and write one network TSV with columns `source`, `target`, `importance`, `pvalue`, `pvalue_bh`. Defaults: `--regulator_types TF SF TF_SF`, `--num_target_clusters 1`, `--num_permutations 1000`. signifikante's kmeans target clustering needs at least 40 targets and 40 samples.

- `whippet_node_grn`: regulators are the informative Whippet nodes of regulator genes (TF, SF, TF_SF) (`informative_exons.informative_nodes`), as Psi named `<gene>:<node>`; targets are the expressed transcripts (TPM). Missing Psi (fewer than `--min_reads`) is filled with the node mean after dropping nodes missing in more than `--max_missing` of the samples. Sample names are the `<sample>.psi(.gz)` basenames and must match the TPM columns.
- `salmon_bootstrap_grn`: transcript → expressed transcript network from Salmon output directories (`salmon_bootstraps.py`, needs `--numBootstraps`). The bootstrap counts are converted to TPM and scaled with the point estimate's StandardScaler, then passed as `bootstraps`. Adds the four SHAP columns.
- `gene_count_grn`: gene → gene network from a merged gene count table (nf-core/rnaseq `salmon.merged.gene_counts.tsv`), normalized to CPM.

The transcript file is a TSV read with `index_col=0`. It has columns `transcript_id`, `gene_id` and one column per sample (TPM).

## Architecture

### Inference ([infer_networks_v3.py](src/alternet/runners/infer_networks_v3.py), [inference.py](src/alternet/inference.py))

Gene expression is the per-`gene_id` sum of transcripts. Both matrices are transposed to samples × features and **standardized** (`standardize_dataframe`) before inference. The three networks:

| File arg | Regulators | Targets | Data |
|---|---|---|---|
| `--canonical` | TF/TF_SF genes | all genes | gene matrix |
| `--source_as` | TF/TF_SF transcripts | all genes | hybrid matrix (`create_hybrid_data`: TF transcript columns + all gene columns) |
| `--fully_as` | all regulator transcripts (TF, SF, TF_SF) | all transcripts | transcript matrix |

`inference()` starts a Dask `LocalCluster` and runs `compute_grn` `n_runs` times. `compute_grn` chunks targets in groups of 5000 and calls `grnboost2`. `aggregate_results` then groups by (source, target) into `frequency`, `mean_importance` and `median_importance`. If you pass `save_dir`, each run is cached as `grn_{i}.tsv`.

### Regulator list

`data/regulator_db.tsv` has columns `Regulator_name`, `Gene stable ID`, `Transcript stable ID`, `Regulator_type` ∈ {`TF`, `SF`, `TF_SF`}. It is built with `annotation.map_tf_ids` / `map_sf_ids` / `combine_tf_sf_lists`. A gene that is in both lists becomes `TF_SF`. The SF list (`data/splicefactors.csv`) merges SpliceAid-F, SplicingLore, KEGG spliceosome and ENCODE RBP sources (see `data/README.md`).

### Downstream ([alternet_pipeline_v3.py](src/alternet/runners/alternet_pipeline_v3.py))

1. `compare_nets.canonical_names` adds `source_type`/`target_type` (`gene` if the ID starts with `ENSG`, otherwise `transcript`), `source_gene`, `target_gene`, `source_transcript`, `target_transcript` and `reg_type`.
2. `compare_nets.filter_edges` keeps edges with `frequency >= 10`, then the top 30% by `median_importance` (the 0.7 quantile), separately per network.
3. The three networks are concatenated. `get_best_variable` groups edges by gene-level key `edge_gg` (`source_gene_target_gene`). It compares importance across edge types (`gene-gene`, `transcript-gene`, `transcript-transcript`) and assigns a `category` of `specific_strict`, `specific`, `equivalent` or `other`. It also sets `is_unique_edge`.
4. `plausibility_filtering` (from `compare_nets`) sets `is_plausible = False` on unique edges when both regulator and target genes are dominated by one isoform (`r_dom=0.9`). It does not drop any rows.
5. The edges are split by regulator type:
   - **TF net**: `reg_type == 'TF'`.
   - **SF net**: SF transcript → transcript edges go to `splicefactor_evidence.compute_set_c`. It tests Spearman correlation between SF expression and target-isoform *usage* (`calculate_transcript_usage`), applies FDR, and requires a minimum delta-usage.
   - **TF_SF**: `tf_sf_disambigouation_fully_as_aware` labels each edge `tfsf_tf_like`, `tfsf_sf_like`, `tfsf_joint` or `tfsf_ambiguous`. `tfsf_tf_like` edges join the TF net, `tfsf_sf_like` edges go through `compute_set_c` into the SF net, and the rest form the ambi net. This step uses joblib with `n_cores=16`.
6. **Annotation**: `annotation.create_transcipt_annotation_database` builds a database from APPRIS + DIGGER. `annotate_isoform_exclusive_edges` runs once for the source transcript and once for the target transcript (suffixes `_source`/`_target`) and attaches the exons/domains each isoform has or lacks compared with its siblings.

### Module status

- Used by v3: `inference.py`, `data_preprocessing.py`, `compare_nets.py`, `splicefactor_evidence.py`, `annotation.py`.
- Not used by any runner: `plausibility_filtering.py` (`filter_set_a`–`filter_set_d`) and `postprocessing.py`. Each defines a function with the same name as one in `compare_nets.py` (`plausibility_filtering`, `filter_edges`), but with a different signature. Make sure you import the right one, since the runners use `from ... import *`.
- `gtex_dataloader.py` holds GTEx loading helpers and is star-imported by `infer_networks_v3`.
- `deprecated/`: old pipeline, tutorials and scripts. Don't add new logic there.

### Column conventions

Identifiers are Ensembl stable IDs without version suffixes. The biomart columns `Gene stable ID` and `Transcript stable ID` are the join keys (constants `GENE_STABLE_ID` / `TRANSCRIPT_STABLE_ID` in `annotation.py`). The mappings `tx2gene` and `gene2tx` (gene → set of transcripts) are built from biomart in the runners.

## Known quirks and bugs

- `infer_networks_v3`: `--n_runs` is parsed but never used (`N_RUNS = 10` is never passed on either). `inference()` runs its default `n_runs=10`.
- The third argument of `canonical_names` is the regulator-type map (`gene2regtype`, keyed by both transcript and gene IDs), not `gene2tx`.
- `alternet_pipeline_v2` calls `filter_edges(..., max_p_value=...)`, which needs `filter_edges_pvalue` instead.
- `data_preprocessing.create_hybrid_data` is defined twice; the second definition wins.
- `compute_dominance_metrics` exists in both `splicefactor_evidence.py` and `plausibility_filtering.py`.
- Misspellings that are part of identifiers: `create_transcipt_annotation_database`, `tf_sf_disambigouation_fully_as_aware`, and `ambi`/`ambigous` in older outputs. Don't rename them in only one place.
- `configs/*.yaml` are path templates for the old pipeline. The v3 runners don't read them.
