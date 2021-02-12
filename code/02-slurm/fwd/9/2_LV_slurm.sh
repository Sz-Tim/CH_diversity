#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24 
#SBATCH --job-name="LV_9-2" 
#SBATCH --output=logs/LV_9-2.output 

module purge 
module load linuxbrew/colsa 

srun bash code/02-slurm/fwd/9/2_LV.sh