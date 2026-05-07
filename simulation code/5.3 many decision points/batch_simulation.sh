#!/bin/bash
#
# batch_simulation.sh
#
# Submit all 100 tasks for Section 5.5 (large-T regime) simulation.
# Grid: 2 sample sizes x 5 T values = 10 settings, x 10 seeds = 100 tasks.
# 1000 MC reps per setting (nsim_per_seed = 100).

n_tasks=100

mkdir -p out_and_err_files

for i in $(seq 1 $n_tasks); do
    sbatch simulation.slurm $i
done
