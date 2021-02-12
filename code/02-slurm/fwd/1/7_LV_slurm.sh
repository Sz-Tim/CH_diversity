#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_1-7" 
#SBATCH --output=logs/LV_1-7.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/1/7_LV.sh