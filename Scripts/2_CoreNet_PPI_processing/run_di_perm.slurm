#!/bin/bash
#SBATCH --job-name=di_perm
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/di_perm_%A_%a.out
#SBATCH --error=logs/di_perm_%A_%a.err
#SBATCH --array=0-4

# ----------------------------------------------------------------------
# Unified DI permutation pipeline. One script, five sets, job array.
#
# Set list (index -> name):
#   0  PPS_noprop      PPS (no propagation) significant genes
#   1  logFC           logFC RNA gene list
#   2  CoreNet         union of CoreNet-PPI and CoreNet-GRN
#   3  CoreNet-GRN     CoreNet-GRN only
#   4  CoreNet-PPI     CoreNet-PPI only
#
# Submit all five:
#     sbatch run_di_perm.slurm
#
# Submit a subset (e.g. only PPS_noprop and CoreNet):
#     sbatch --array=0,2 run_di_perm.slurm
#
# Submit just one (e.g. CoreNet-GRN):
#     sbatch --array=3 run_di_perm.slurm
# ----------------------------------------------------------------------

set -euo pipefail
mkdir -p logs

echo "JobID:        ${SLURM_JOB_ID}"
echo "ArrayJobID:   ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "ArrayTaskID:  ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Node:         $(hostname)"
echo "Started:      $(date)"
echo "Working dir:  $(pwd)"
echo "------------------------------------------------------------"

# --- EDIT-ME: shared paths -------------------------------------------
SCRIPT="first_degree_interactors_pipeline.py"

HINT_H5="/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5"
ALL_HINT_GENES="../../Data/all_hint_genes_list.pkl"

# Literature sets are now CSVs (single 'Genes' column).
VINUESA="../../Data/literature_sets/vinuesa_list.csv"
HART="../../Data/literature_sets/TFH_review_list.csv"
CROTTY="../../Data/literature_sets/crotty_list.csv"

# Gold-standard is also a CSV now (single 'Genes' column).
GOLD_STANDARD="../../Data/literature_sets/vinuesa_list.csv"

# Top-level output dir; each set writes to ${BASE_OUTPUT_DIR}/${SET_NAME}/
BASE_OUTPUT_DIR="../Figure2/outputs/upset_plot_overlaps"

# --- shared run params ------------------------------------------------
NUM_PERM=1000
MODE="all"                       # 'di_then_match' | 'match_then_di' | 'sig_only' | 'both' | 'all' | 'none'
SAVE_FORMATS="csv"
SEED=42

TRUE_DI=1
SAVE_CYTOSCAPE=1
CYTOSCAPE_VARIANT="cytoscape"    # 'simple' | 'goldstand' | 'cytoscape' | 'all'

# --- per-set definitions ----------------------------------------------
# Index in each array corresponds to the SLURM_ARRAY_TASK_ID.
SET_NAMES=(
    "PPS_noprop"
    "logFC"
    "CoreNet"
    "CoreNet_GRN"
    "CoreNet_PPI"
)

# Input files. For CoreNet (index 2) this is built on the fly from
# the two component files below; the entry here is a placeholder.
SIG_GENES_FILES=(
    "../../Data/EBSeq_Genes_top184.csv"
    "../../Data/EBSeq_RNASeq_logFC_top184_0.01.csv"
    "__BUILD_CORENET_UNION__"
    "../../Data/GRN_genes_all_sets_152.csv"
    "../../Data/PPS_1_significant_gene_list.csv"
)

# Column name holding gene symbols in each input CSV.
GENE_COLS=(
    "Gene"
    "gene"
    "Genes"
    "Genes"
    "Genes"
)

# Prefix used in output filenames (Python doesn't love hyphens in filenames
# from a downstream-tooling perspective, but they're fine for the pipeline).
PREFIXES=(
    "PPS_noprop_184"
    "logFC"
    "CoreNet"
    "CoreNet_GRN"
    "CoreNet_PPI"
)

# CoreNet (index 2) is built as the union of these two:
CORENET_PPI_INPUT="../../Data/PPS_1_significant_gene_list.csv"
CORENET_GRN_INPUT="../../Data/GRN_genes_all_sets_152.csv"

# --- pick the set for this array task ---------------------------------
IDX="${SLURM_ARRAY_TASK_ID:-0}"

if [[ "${IDX}" -ge "${#SET_NAMES[@]}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${IDX} is out of range (max $((${#SET_NAMES[@]} - 1)))"
    exit 1
fi

SET_NAME="${SET_NAMES[${IDX}]}"
SIG_GENES="${SIG_GENES_FILES[${IDX}]}"
GENE_COL="${GENE_COLS[${IDX}]}"
PREFIX="${PREFIXES[${IDX}]}"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${SET_NAME}"

mkdir -p "${OUTPUT_DIR}"

echo "Set:          ${SET_NAME}  (array idx ${IDX})"
echo "Prefix:       ${PREFIX}"
echo "Output dir:   ${OUTPUT_DIR}"
echo "------------------------------------------------------------"

# --- build the CoreNet union on the fly if needed ---------------------
if [[ "${SIG_GENES}" == "__BUILD_CORENET_UNION__" ]]; then
    SIG_GENES="${OUTPUT_DIR}/CoreNet_union_genes.csv"
    echo "Building CoreNet union from CoreNet_PPI + CoreNet_GRN..."
    python - <<PYEOF
import pandas as pd
ppi = pd.read_csv("${CORENET_PPI_INPUT}")
grn = pd.read_csv("${CORENET_GRN_INPUT}")
union = sorted(set(ppi['Genes']) | set(grn['Genes']))
pd.DataFrame({'Genes': union}).to_csv("${SIG_GENES}", index=False)
print(f"Union size: {len(union)} (PPI={len(set(ppi['Genes']))}, GRN={len(set(grn['Genes']))})")
PYEOF
    echo "------------------------------------------------------------"
fi

# --- assemble args ----------------------------------------------------
ARGS=(
    --sig_genes      "${SIG_GENES}"
    --gene_col       "${GENE_COL}"
    --prefix         "${PREFIX}"
    --hint_h5        "${HINT_H5}"
    --all_hint_genes "${ALL_HINT_GENES}"
    --vinuesa        "${VINUESA}"
    --hart           "${HART}"
    --crotty         "${CROTTY}"
    --output_dir     "${OUTPUT_DIR}"
    --num_perm       "${NUM_PERM}"
    --mode           "${MODE}"
    --save_formats   "${SAVE_FORMATS}"
    --seed           "${SEED}"
)

[[ "${TRUE_DI}" == "1" ]] && ARGS+=(--save_di)

if [[ "${SAVE_CYTOSCAPE}" == "1" ]]; then
    ARGS+=(--save_cytoscape --cytoscape_variant "${CYTOSCAPE_VARIANT}")
    [[ -n "${GOLD_STANDARD}" ]] && ARGS+=(--gold_standard "${GOLD_STANDARD}")
fi

echo "Command:"
echo "  python ${SCRIPT} ${ARGS[*]}"
echo "------------------------------------------------------------"

python "${SCRIPT}" "${ARGS[@]}"

echo "------------------------------------------------------------"
echo "Finished:     $(date)"
