#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_2-5" 
#SBATCH --output=logs/LV_2-5.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/2/5_LV.sh