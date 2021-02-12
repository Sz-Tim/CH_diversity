#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_10-2" 
#SBATCH --output=logs/cov_10-2.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/10/2_cov.sh