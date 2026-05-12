# Wall-clock benchmark

Reproduces Table F.1 in the revised paper (Appendix F.1, "Practical
Implementation of the Proposed Algorithms"). Times Algorithm 1 of the main
paper (no cross-fitting) and Algorithm D.1 of the appendix (cross-fitting,
K = 5) for the continuous proximal outcome on the Drink Less MRT dataset
(n = 349, T = 30), single-threaded on a local machine. Added in response to
Reviewer 1 comment 3 on computational scalability.

## Files

- `run_benchmark.R` — main script. Times the four nuisance-fit methods (WCLS,
  GAM, RF, super learner) with and without cross-fitting.
- `functions/` — copies of the estimator code (`wcls_eif.R`, `eif.R`,
  `fit_d.R`, `wcls_original.R`).
- `result/` — output from a prior run on the manuscript author's machine
  (`wallclock_continuous.csv`, `wallclock_continuous.RDS`, `run_log.txt`).

## Data

`run_benchmark.R` expects `data/DrinkLess_cont_bin.RDS` (a long-format data
frame with columns `ID`, `dp`, `A`, `prob_A`, `Y_cont`, etc.). The continuous
outcome version of the Drink Less data is **not publicly available**; see the
top-level `README.md` of this repository for guidance on data access.

## How to run

```r
# In RStudio, or set the working directory to this folder
source("run_benchmark.R")
```
