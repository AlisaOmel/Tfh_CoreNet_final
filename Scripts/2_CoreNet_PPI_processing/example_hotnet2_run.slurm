#!/bin/bash

#SBATCH -t 6-00:00 
#SBATCH --job-name=B_hotnet_TFH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4

module purge
module load python/anaconda2.7-4.2.0

# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
#hotnet2=..
num_cores=-1

python ./HotNet2.py \
-nf  ../../Input_data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.h5 \
-pnp ../../Input_data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4_##NUM##.h5 \
-hf  ../../Data/EBSeq_Genes_Scores.txt \
-np 500 \
-hp 500 \
-o ../../Sample_outputs/PPI_network \
-c  $num_cores
