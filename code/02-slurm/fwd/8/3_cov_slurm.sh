#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_8-3" 
#SBATCH --output=logs/cov_8-3.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/8/3_cov.sh