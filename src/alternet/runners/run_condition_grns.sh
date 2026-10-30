#!/usr/bin/env bash
# Splits a gene expression folder by condition and runs the three single-network wrappers
# (whippet_node_grn, gene_count_grn, salmon_bootstrap_grn) on every condition.
#
# Usage: bash run_condition_grns.sh <data_dir> <mapping.tsv> <output_dir>
#   data_dir:   gene expression folder (<SRR>/counts/{salmon,whippet}/, merged tables in counts/)
#   mapping:    TSV with a header, sample ID and condition (e.g. data/GSE278521_mapping.tsv)
#   output_dir: gets <condition>/ (moved sample directories, split tables) and networks/
#
# Set the reference files and options through these environment variables:
#   ALTERNET_DIR   alternet repository (default: two levels above this script's package)
#   GTF            GTF used for the Whippet index
#   EXON_NODES     Whippet graph.exons.tab(.gz)
#   PYTHON         Python command with alternet installed (default: python)
#   EXTRA_ARGS     extra arguments for all three wrappers (e.g. "--num_permutations 100")
# A network is computed only if its output file doesn't exist yet. The split is skipped if
# output_dir/samples.tsv exists (the sample directories have been moved already).
set -euo pipefail

DATA_DIR=$1
MAPPING=$2
OUT_DIR=$3

ALTERNET_DIR=${ALTERNET_DIR:-$(cd "$(dirname "$0")/../../.." && pwd)}
GTF=${GTF:?Set GTF to the GTF used for the Whippet index}
EXON_NODES=${EXON_NODES:?Set EXON_NODES to the Whippet graph.exons.tab(.gz)}
PYTHON=${PYTHON:-python}
EXTRA_ARGS=${EXTRA_ARGS:-}
REGULATORS=$ALTERNET_DIR/data/regulator_db.tsv
BIOMART=$ALTERNET_DIR/data/biomart.txt

if [[ ! -f $OUT_DIR/samples.tsv ]]; then
    $PYTHON -m alternet.runners.split_by_condition --data_dir "$DATA_DIR" --mapping "$MAPPING" --output_dir "$OUT_DIR"
fi

NET_DIR=$OUT_DIR/networks
mkdir -p "$NET_DIR"

for COND_DIR in "$OUT_DIR"/*/; do
    COND=$(basename "$COND_DIR")
    [[ $COND == networks ]] && continue
    TPM_FILE=$(ls "$COND_DIR"counts/*_salmon_transcript_tpm.csv)
    echo "=== $COND: $(ls -d "$COND_DIR"*/counts | wc -l) samples"

    OUT=$NET_DIR/$COND.whippet_node_grn.tsv
    if [[ ! -f $OUT ]]; then
        $PYTHON -m alternet.runners.whippet_node_grn \
            --psi_files "$COND_DIR"*/counts/whippet/*.psi.gz --tpm_file "$TPM_FILE" \
            --gtf "$GTF" --exon_nodes "$EXON_NODES" --regulator_file "$REGULATORS" \
            --output "$OUT" --nodes_output "$NET_DIR/$COND.whippet_nodes.tsv" \
            --benchmark "$NET_DIR/$COND.whippet_node_grn.benchmark.tsv" $EXTRA_ARGS
    fi

    OUT=$NET_DIR/$COND.gene_count_grn.tsv
    if [[ ! -f $OUT ]]; then
        $PYTHON -m alternet.runners.gene_count_grn \
            --tpm_file "$TPM_FILE" --biomart_file "$BIOMART" --regulator_file "$REGULATORS" \
            --output "$OUT" --benchmark "$NET_DIR/$COND.gene_count_grn.benchmark.tsv" $EXTRA_ARGS
    fi

    OUT=$NET_DIR/$COND.salmon_bootstrap_grn.tsv
    if [[ ! -f $OUT ]]; then
        N_QUANT=$(find "$COND_DIR" -mindepth 3 -maxdepth 3 -path '*/counts/salmon' -type d | wc -l)
        N_BOOT=$(find "$COND_DIR" -mindepth 6 -maxdepth 6 -path '*/counts/salmon/aux_info/bootstrap/bootstraps.gz' | wc -l)
        if [[ $N_BOOT -eq $N_QUANT ]]; then
            $PYTHON -m alternet.runners.salmon_bootstrap_grn \
                --quant_dirs "$COND_DIR"*/counts/salmon --regulator_file "$REGULATORS" \
                --output "$OUT" --benchmark "$NET_DIR/$COND.salmon_bootstrap_grn.benchmark.tsv" $EXTRA_ARGS
        else
            echo "Skipping salmon_bootstrap_grn for $COND: $N_BOOT of $N_QUANT samples have Salmon bootstraps (salmon quant --numBootstraps)"
        fi
    fi
done
