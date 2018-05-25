#!/bin/bash


#SBATCH --job-name=THIS_EXPERIMENT
#SBATCH --output=THIS_EXPERIMENT_OUTPUT.txt
#SBATCH --error=THIS_EXPERIMENT_ERROR.txt 
#SBATCH --time=36:00:00
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000


Rscript THIS_EXPERIMENT_R_SCRIPT.R $SLURM_ARRAY_TASK_ID 
