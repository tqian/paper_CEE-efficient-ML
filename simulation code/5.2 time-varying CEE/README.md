# Simulation 5.4 — Time-Varying CEE (EJS revision)

This directory contains the code for the new Section 5.4 of the manuscript,
added in response to Reviewer 2's minor comment 2: *"In simulations, it will be
interesting to see some settings where $\gamma_t$ depends on $t$."*

## Generative model

Continuous outcome, $T = 10$. Same skeleton as the original
`dgm_cont_time_s_tvvar.R` from `original submission - github code/`, with one
change: the CEE is now time-moderated,
$$
\gamma_t(\beta) = \beta_0 + \beta_1 \cdot (t / T),
$$
with $\beta_0 = 0.5$, $\beta_1 = -0.3$ (effect $\approx 0.47$ at $t=1$, $0.20$
at $t=10$). The covariate $S_t \sim \mathrm{Unif}[-2, 2]$ is retained as a
driver of $\mu_t(H_t, 0)$ via the unchanged `linear / sine / dbeta / step`
patterns.

## Estimators

Four estimators run inside each Monte Carlo iteration on the same generated
dataset (paired Monte Carlo for tighter relative-efficiency comparisons):

1. **WCLS** — baseline, no augmentation
2. **GAM (cross-fitted)** — $\hat\mu_t$ fit by `mgcv::gam`, 10 folds
3. **SL.smooth (cross-fitted)** — $\hat\mu_t$ fit by `SuperLearner` with
   smooth library (`SL.mean, SL.glm, SL.gam, SL.earth`), 10 folds
4. **Oracle** — true $\mu_t$ (from the dgm), empirical $d_t$

In all four, the moderator is `dp_norm = t/T` and the working covariates for
$\hat\mu_t$ are `c("dp", "S")`.

## Sweep grid

Mirrors Figure 1 of the existing manuscript:

- left panel:   $n \in \{30, 50, 100\}$ at $\lambda_1 = 1$
- middle/right: $\lambda_1 \in \{0, 0.5, 1, 1.5, 2, 2.5, 3\}$ at $n = 100$
- 4 control patterns: `linear`, `sine`, `dbeta`, `step`

Total: 36 unique settings, 1000 MC reps each, packed as
`nsim_per_seed = 100 x n_seeds = 10` $\Rightarrow$ **360 SLURM tasks**.

## Workflow

```bash
# 1. submit all 360 tasks
bash batch_simulation.sh

# 2. (optional) check progress
bash check_missing.sh | wc -l

# 3. resubmit any failed tasks
bash batch_rerun_missing.sh

# 4. once all tasks finish, aggregate results
sbatch collect_results.slurm    # or: Rscript collect_results.R

# 5. build figures
sbatch make_plot.slurm          # or: Rscript make_plot.R
```

## File layout

```
function/                       # helpers (copied unchanged from original sim except the dgm)
    dgm_cont_time_varying_CEE.R   # NEW: time-moderated CEE
    wcls_eif.R                    # cross-fitted EIF estimator
    eif.R                         # core EIF root-finder + variance
    fit_d.R                       # data-adaptive d_t
    wcls_original.R               # WCLS baseline

simulation_R_code.R           # main driver, takes <itask> as arg
simulation.slurm              # SLURM wrapper for one task
batch_simulation.sh           # submits all 360 tasks
check_missing.sh              # lists tasks with no result_tmp/<itask>.RDS
batch_rerun_missing.sh        # resubmits only the missing ones

collect_results.R             # aggregates result_tmp/ into result_collected/
collect_results.slurm

make_plot.R                   # builds Section 5.4 figure
make_plot.slurm

result_tmp/                   # per-task RDS files, populated by simulation
result_collected/             # tidy aggregated results, populated by collect_results.R
plot/                         # figures, populated by make_plot.R
out_and_err_files/            # SLURM stdout/stderr
```
