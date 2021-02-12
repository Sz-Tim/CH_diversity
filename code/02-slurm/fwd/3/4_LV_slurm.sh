#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_3-4" 
#SBATCH --output=logs/LV_3-4.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/3/4_LV.sh