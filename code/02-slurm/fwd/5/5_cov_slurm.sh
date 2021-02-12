#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_5-5" 
#SBATCH --output=logs/cov_5-5.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/5/5_cov.sh