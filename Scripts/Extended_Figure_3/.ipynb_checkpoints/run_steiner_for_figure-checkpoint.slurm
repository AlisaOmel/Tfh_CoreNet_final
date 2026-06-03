#!/bin/bash
#SBATCH --job-name=run_steiner
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster=smp
#SBATCH --mem 265G

python steiner_tree_extfig.py

echo done
