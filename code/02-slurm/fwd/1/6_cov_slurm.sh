#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_1-6" 
#SBATCH --output=logs/cov_1-6.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/1/6_cov.sh