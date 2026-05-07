# Simulation 5.5 — Large-T regime (EJS revision)

This directory contains the code for the new Section 5.5 of the manuscript,
added in response to Reviewer 2's minor comment 2: *"...when [$\gamma_t$ does
not depend on $t$], it will be interesting to see results when $n$ is fixed
and $T$ is large."*

## Generative model

Continuous outcome, **constant CEE** $\gamma_t(\beta) = \beta_0 = 0.5$.
Same skeleton as `dgm_cont_time_s_tvvar.R` from `original submission - github
code/`, with the marginal CEE form. $S_t \sim \mathrm{Unif}[-2, 2]$ is retained
as a covariate driving $\mu_t(H_t, 0)$ via the unchanged `dbeta` pattern (a
smooth unimodal bump on $t/T \in [0,1]$, whose shape does not get harder to
estimate as $T$ grows).

## Estimators

Same four as sim 5.4 (paired Monte Carlo per iteration):

1. **WCLS** — baseline
2. **GAM (cross-fitted)**
3. **SL.smooth (cross-fitted)**
4. **Oracle**

In all four, `moderator = NULL` (marginal CEE) and `control = c("dp", "S")`.

## Sweep grid

- $n \in \{30, 100\}$
- $T \in \{10, 30, 50, 100, 200\}$
- `control_pattern = "dbeta"` with $\lambda_1 = 1$ (fixed)
- $\rho = 0.5$

Total: 10 settings, 1000 MC reps each, packed as
`nsim_per_seed = 100 x n_seeds = 10` $\Rightarrow$ **100 SLURM tasks**.

## Workflow

```bash
bash batch_simulation.sh             # submit 100 tasks
bash check_missing.sh | wc -l        # progress
bash batch_rerun_missing.sh          # if any failed
sbatch collect_results.slurm
sbatch make_plot.slurm
```
