#!/bin/bash

#SBATCH --job-name=sim_long_model
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt 
#SBATCH --array=1-150
#SBATCH --time=10:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=MaxMemPerCPU

Rscript fit_probit_model.R $SLURM_ARRAY_TASK_ID