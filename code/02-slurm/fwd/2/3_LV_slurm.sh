#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_2-3" 
#SBATCH --output=logs/LV_2-3.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/2/3_LV.sh