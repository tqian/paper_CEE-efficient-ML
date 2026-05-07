#!/bin/bash
#
# check_missing.sh
#
# Print task IDs (1..100) whose result_tmp/<itask>.RDS file is missing.

n_tasks=100

for i in $(seq 1 $n_tasks); do
    if [ ! -f "result_tmp/${i}.RDS" ]; then
        echo $i
    fi
done
