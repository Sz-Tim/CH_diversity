#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_1-2" 
#SBATCH --output=logs/LV_1-2.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/1/2_LV.sh