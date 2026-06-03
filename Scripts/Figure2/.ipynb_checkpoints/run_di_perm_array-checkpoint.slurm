#!/bin/bash
#SBATCH --job-name=di_perm
#SBATCH --partition=htc
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/di_perm_%A_%a.out
#SBATCH --error=logs/di_perm_%A_%a.err

# ----------------------------------------------------------------------
# Direct-interactor permutation pipeline -- SLURM array over 4 gene sets.
#
# Tasks (SLURM_ARRAY_TASK_ID):
#   0 = CoreNet_PPI       (from PPS_1_significant_gene_list.csv)
#   1 = CoreNet_GRN       (from final_genes_all_sets_152.csv)
#   2 = CoreNet_PPI_GRN   (union of the two above; built inline)
#   3 = logFC_rna         (from Pattern_1234_RNASeq_S1234_logFC_top184_0.01_correct.csv)
#
# Submit:  sbatch run_di_perm_array.slurm
# All four jobs run concurrently subject to QoS limits.
# ----------------------------------------------------------------------

set -euo pipefail
mkdir -p logs

echo "JobID:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node:         $(hostname)"
echo "Started:      $(date)"
echo "Working dir:  $(pwd)"
echo "------------------------------------------------------------"


# --- shared paths -----------------------------------------------------
SCRIPT="first_degree_interactors_pipeline.py" 
CORENET_PPI_INPUT="/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/PPS_1_significant_gene_list.csv"
CORENET_GRN_INPUT="/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/processed_Taiji/final_genes_all_sets_152.csv"
LOGFC_INPUT="/ix3/djishnu/Alisa/Tfh/no_propagation/Pattern_1234_RNASeq_S1234_logFC_top184_0.01_correct.csv"

HINT_H5="/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5"
ALL_HINT_GENES="/ix/djishnu/Alisa/Tfh/Network_analysis/results/all_hint_genes_list.pkl"

VINUESA="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/vinuesa_list.pkl"
HART="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/TFH_review_list.pkl"
CROTTY="/ix/djishnu/Alisa/Tfh/ForPaper/literature_sets/mouse_gene_list.pkl"



# --- run params (shared across all tasks) -----------------------------
NUM_PERM=1000
MODE="all"              # 'di_then_match' | 'match_then_di' | 'sig_only' | 'both' | 'all' | 'none'
SAVE_FORMATS="csv"
SEED=42
TRUE_DI=1
SAVE_CYTOSCAPE=0
CYTOSCAPE_VARIANT="simple"
GOLD_STANDARD=""

# --- per-task config --------------------------------------------------
# Each task picks its sig_genes file, gene_col, and prefix.
# The CoreNet_PPI_GRN task computes the union inline below before running.
case "${SLURM_ARRAY_TASK_ID}" in
    0)
        PREFIX="CoreNet_PPI"
        SIG_GENES="${CORENET_PPI_INPUT}"
        GENE_COL="Genes"
        ;;
    1)
        PREFIX="CoreNet_GRN"
        SIG_GENES="${CORENET_GRN_INPUT}"
        GENE_COL="Genes"
        ;;
    2)
        PREFIX="CoreNet_PPI_GRN"
        # Build the union inline. Write to a per-task scratch dir.
        UNION_DIR="${OUTPUT_DIR}/_union_inputs"
        mkdir -p "${UNION_DIR}"
        SIG_GENES="${UNION_DIR}/CoreNet_PPI_GRN_union.csv"
        GENE_COL="Genes"

        echo "Building CoreNet_PPI_GRN union from CoreNet_PPI + CoreNet_GRN..."
        python - <<PYEOF
import pandas as pd
ppi = pd.read_csv("${CORENET_PPI_INPUT}")
grn = pd.read_csv("${CORENET_GRN_INPUT}")
union = sorted(set(ppi['Genes']) | set(grn['Genes']))
pd.DataFrame({'Genes': union}).to_csv("${SIG_GENES}", index=False)
print(f"Union size: {len(union)} (PPI={len(set(ppi['Genes']))}, GRN={len(set(grn['Genes']))})")
PYEOF
        ;;
    3)
        PREFIX="logFC_rna"
        SIG_GENES="${LOGFC_INPUT}"
        GENE_COL="gene"
        ;;
    *)
        echo "Unknown SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}" >&2
        exit 1
        ;;
esac

echo "Task ${SLURM_ARRAY_TASK_ID}: prefix=${PREFIX}, sig_genes=${SIG_GENES}, gene_col=${GENE_COL}"
echo "------------------------------------------------------------"

OUTPUT_DIR="/ix/djishnu/Alisa/Tfh/testing_for_github/Figure2/test_output/${PREFIX}/"
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
