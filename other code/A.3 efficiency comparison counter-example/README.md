# Efficiency comparison counter-example (Appendix A.3)

Empirical Monte Carlo verification of the counter-example in Appendix A.3 of
the revised paper, where carryover from past treatment makes the proposed
covariate-adjusted estimators *less* efficient than the unadjusted IPW
estimator.

Added in response to Reviewer 2 comment 3 ("Theoretical guarantees"). The
sharper version (`sim_counterexample_v2.R`) is the one cited in the paper.

## Files

- `sim_counterexample.R` — initial exploration. T = 2 with carryover via
  `Y_3 = alpha * (Y_2 - 0.5) * (A_2 - 0.5) + A_2 + eps_2`. Shows that
  adjustment with the true conditional mean inflates the asymptotic variance
  of the unweighted (constant-d) version.
- `sim_counterexample_v2.R` — sharper counter-example. T = 2 with carryover
  via `Y_3 = theta * (A_1 - 0.5) * (A_2 - 0.5) + A_2 + eps_2`, theta = -3.
  Shows that even the *optimal* proposed estimator (with d_t^* and the true
  mu_t^*) is asymptotically less efficient than the unadjusted estimator.
  This is the version cited in Section 3.3 / Appendix A.3 of the manuscript.

## How to run

```r
# In RStudio, or set the working directory to this folder
source("sim_counterexample_v2.R")
```

`sim_counterexample_v2.R` sources `eif.R` and `wcls_original.R` from
`../D.1 wallclock benchmark/functions/`. `sim_counterexample.R` is
self-contained.
