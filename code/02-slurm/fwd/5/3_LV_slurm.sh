#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_5-3" 
#SBATCH --output=logs/LV_5-3.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/5/3_LV.sh