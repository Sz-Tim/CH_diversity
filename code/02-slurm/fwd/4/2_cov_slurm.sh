#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_4-2" 
#SBATCH --output=logs/cov_4-2.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/4/2_cov.sh