#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="cov_2-7" 
#SBATCH --output=logs/cov_2-7.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/2/7_cov.sh