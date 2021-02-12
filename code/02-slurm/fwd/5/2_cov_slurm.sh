#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_5-2" 
#SBATCH --output=logs/cov_5-2.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/5/2_cov.sh