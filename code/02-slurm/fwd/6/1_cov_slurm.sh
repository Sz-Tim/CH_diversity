#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_6-1" 
#SBATCH --output=logs/cov_6-1.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/6/1_cov.sh