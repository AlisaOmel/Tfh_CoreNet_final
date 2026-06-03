#!/bin/bash
#SBATCH --job-name=RNA_ATAC
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process
#SBATCH --cluster=smp

module purge
module load taiji/1.2.0

cd ATAC_RNA

export PATH=$PATH:/softwares/bwa-0.7.17
export PATH=$PATH:/softwares/samtools-1.16.1

taiji run --config config_atac.yml

echo done
