#!/bin/bash

# batch_simulation.sh
#
# Submit all 360 tasks for Section 5.4 (time-varying CEE) simulation.
# Grid:
#   9 (sample_size, lambda_1) combos x 4 control_patterns = 36 settings
#   x 10 seeds = 360 tasks
#   1000 MC reps per setting (nsim_per_seed = 100)

n_tasks=360

mkdir -p out_and_err_files

for i in $(seq 1 $n_tasks); do
    sbatch simulation.slurm $i
done
