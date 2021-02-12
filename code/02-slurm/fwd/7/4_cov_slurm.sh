#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_7-4" 
#SBATCH --output=logs/cov_7-4.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/7/4_cov.sh