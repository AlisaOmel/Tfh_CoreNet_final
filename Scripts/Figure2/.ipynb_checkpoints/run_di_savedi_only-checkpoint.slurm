#!/bin/bash
#SBATCH --job-name=di_savedi_only
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/di_savedi_only_%j.out
#SBATCH --error=logs/di_savedi_only_%j.err

# ----------------------------------------------------------------------
# Run ONLY the --save_di step for PPS_noprop_184.
# No permutations -- this just computes and writes the real DI gene set.
# Should finish in well under a minute.
# ----------------------------------------------------------------------

set -euo pipefail

mkdir -p logs

echo "JobID:        ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "Started:      $(date)"
echo "Working dir:  $(pwd)"
echo "------------------------------------------------------------"


# --- EDIT-ME: paths ---------------------------------------------------
SCRIPT="/path/to/first_degree_interactors_pipeline.py"

SIG_GENES="/ix3/djishnu/Alisa/Tfh/no_propagation/Pattern_1234_Genes_top184_correct.csv"
GENE_COL="Gene"
PREFIX="PPS_noprop_184"

HINT_H5="/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5"
ALL_HINT_GENES="/ix/djishnu/Alisa/Tfh/Network_analysis/results/all_hint_genes_list.pkl"

VINUESA="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/vinuesa_list.pkl"
HART="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/TFH_review_list.pkl"
CROTTY="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/mouse_gene_list.pkl"

# Match the per-prefix folder structure the R script expects
OUTPUT_DIR="/ix/djishnu/Alisa/Tfh/testing_for_github/Figure2/test_output/PPS_noprop_184"

# --- run params -------------------------------------------------------
# MODE='none' skips every permutation pipeline; --save_di is the only work done.
NUM_PERM=1000            # ignored since MODE=none, but argparse requires it
MODE="none"
SAVE_FORMATS="csv"
SEED=42

# Toggles: only TRUE_DI is on
TRUE_DI=1
SAVE_CYTOSCAPE=0
CYTOSCAPE_VARIANT="simple"
GOLD_STANDARD=""

# ----------------------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"

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
