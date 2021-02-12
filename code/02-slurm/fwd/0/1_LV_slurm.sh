#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_0-1" 
#SBATCH --output=logs/LV_0-1.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/0/1_LV.sh