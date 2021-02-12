#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_1-3" 
#SBATCH --output=logs/cov_1-3.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/1/3_cov.sh