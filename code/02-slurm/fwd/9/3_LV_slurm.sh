#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_9-3" 
#SBATCH --output=logs/LV_9-3.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/9/3_LV.sh