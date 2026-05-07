#!/bin/bash
#
# batch_rerun_missing.sh
#
# Resubmit only the tasks whose result_tmp/<itask>.RDS files are missing.
# Usage:
#   bash batch_rerun_missing.sh

missing=$(bash check_missing.sh)
n_missing=$(echo "$missing" | grep -c .)

if [ "$n_missing" -eq 0 ]; then
    echo "All result_tmp files present. Nothing to rerun."
    exit 0
fi

echo "Found $n_missing missing tasks. Resubmitting..."

mkdir -p out_and_err_files

for i in $missing; do
    echo "  Submitting task $i"
    sbatch simulation.slurm $i
done

echo "Done. Submitted $n_missing jobs."
echo "After all jobs complete, run: Rscript collect_results.R"
