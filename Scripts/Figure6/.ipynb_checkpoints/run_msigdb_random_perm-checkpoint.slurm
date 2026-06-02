#!/bin/bash
#SBATCH --job-name=pid_msigdb_perm
#SBATCH -N 1                        # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                  # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1           # Request that ncpus be allocated per process
#SBATCH --output=logs/%x_%j.out     # %x = job name, %j = job ID

# -----------------------------------------------------------------------------
# Usage:
#   sbatch run_permutations.slurm hallmark
#   sbatch run_permutations.slurm pid
#
# Defaults to "hallmark" if no argument is given.
# -----------------------------------------------------------------------------


mkdir -p logs

echo "Running database"
Rscript msigdb_random_permutations_fig6a.R pid
echo done
