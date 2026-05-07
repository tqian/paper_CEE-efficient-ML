#!/bin/bash

#SBATCH --job-name=sim_5_5
#SBATCH -p standard
#SBATCH -A tiancq1_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --error=out_and_err_files/%x.%A.err
#SBATCH --out=out_and_err_files/%x.%A.out

module load R/4.3.3
# $1 = itask (required); $2 = nsim_per_seed (optional, default 100)
Rscript simulation_R_code.R $1 $2
