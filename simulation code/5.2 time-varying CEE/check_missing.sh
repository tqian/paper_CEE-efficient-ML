#!/bin/bash
#
# check_missing.sh
#
# Print task IDs (1..360) whose result_tmp/<itask>.RDS file is missing.
# Usage:
#   bash check_missing.sh            # print missing task IDs
#   bash check_missing.sh | wc -l    # count missing tasks

n_tasks=360

for i in $(seq 1 $n_tasks); do
    if [ ! -f "result_tmp/${i}.RDS" ]; then
        echo $i
    fi
done
