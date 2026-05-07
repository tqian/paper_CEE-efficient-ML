#!/bin/bash

#SBATCH --job-name=sim_5_4         ## job name
#SBATCH -p standard                 ## HPC partition to use
#SBATCH -A tiancq1_lab              ## charge to tiancq1's lab account
#SBATCH --nodes=1                   ## use 1 node, don't ask for multiple
#SBATCH --ntasks=1                  ## number of tasks to launch
#SBATCH --cpus-per-task=1           ## number of cores the job needs
#SBATCH --mem-per-cpu=4G            ## 4 GB per CPU
#SBATCH --time=10:00:00             ## per-task wall time cap
#SBATCH --error=out_and_err_files/%x.%A.err
#SBATCH --out=out_and_err_files/%x.%A.out

module load R/4.3.3
# $1 = itask (required); $2 = nsim_per_seed (optional, default 100)
Rscript simulation_R_code.R $1 $2
