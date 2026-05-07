#!/bin/bash

#SBATCH --job-name=collect_5_4
#SBATCH -p standard
#SBATCH -A tiancq1_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --error=out_and_err_files/%x.%A.err
#SBATCH --out=out_and_err_files/%x.%A.out

module load R/4.3.3
Rscript collect_results.R
