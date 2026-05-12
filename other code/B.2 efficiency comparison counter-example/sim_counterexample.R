# =============================================================================
# Empirical verification of the counter-example in r2c3-theorem-draft.tex,
# Example 6.1 (T = 2 with carryover effect).
#
# DGP:
#   T = 2, simple randomization p_t = 0.5, I_t = 1.
#   Y_2 = A_1 + eps_1
#   Y_3 = alpha * (Y_2 - 0.5) * (A_2 - 0.5) + A_2 + eps_2
#   Working CEE: gamma_t(emptyset; beta) = beta. True beta_star = 1.
#   For alpha < -2, B(0) - B(mu_star) = 2 + alpha < 0, so adjustment with
#   the true conditional mean inflates asymptotic variance for the
#   unweighted (d_t = constant) version.
#
# Estimators compared:
#   A. Unadjusted          : d_t = 1, mu_t = 0          (no covariate adjustment)
#   B. WCLS (Boruvka 2018) : intercept-only OLS on stacked (i, t) data
#   C. DR-WCLS             : d_t = 1/4, mu_t = mu_t^*   (proposed with constant d)
#   D. Optimal proposed    : d_t = d_t^*, mu_t = mu_t^* (proposed with optimal d)
#
# Important construction-specific fact:
#   mu_2^*(H_2, 1) = alpha (Y_2 - 0.5) * 0.5 + 1
#   mu_2^*(H_2, 0) = -alpha (Y_2 - 0.5) * 0.5
# average to b_2^* = 0.5 (constant), and b_1^* = 0.5 too.  Hence the
# WCLS-implied intercept alpha^* = E[Y_{t+1}] = 0.5 coincides with b_t^*,
# so WCLS and DR-WCLS have IDENTICAL asymptotic variance in this construction.
#
# What to expect (theory, with alpha = -3, sigma^2 = 1):
#   V_A (Unadjusted)        = 18.25 / 4 = 4.5625
#   V_B (WCLS)              = 19.25 / 4 = 4.8125  (= V_C)
#   V_C (DR-WCLS)           = 19.25 / 4 = 4.8125
#   V_D (Optimal proposed)  = 1 / (1/4 + 1/15.25) = 3.169
#
# So unadjusted beats WCLS and DR-WCLS by ~5%, but the optimal-d proposed
# estimator beats unadjusted by ~30%.  The counter-example shows that
# fixed-d adjustment can hurt; it does NOT show that the full proposed
# (with optimal d) hurts.
# =============================================================================

set.seed(20260501)

# ----- Configuration ---------------------------------------------------------

alpha   <- -3
sigma   <- 1
n_obs   <- 2000
n_reps  <- 5000
beta_star <- 1

# ----- Theoretical reference values ------------------------------------------

Var_phi1_mu <- 4 * sigma^2
Var_phi2_mu <- alpha^2 * (0.25 + sigma^2) + 4 * sigma^2

V_unadj  <- (2 + 8*sigma^2 + alpha^2*(0.25+sigma^2) + alpha) / 4   # B(0) / 4
V_drwcls <- (8*sigma^2 + alpha^2*(0.25+sigma^2)) / 4               # B(mu*) / 4
V_opt    <- 1 / (1/Var_phi1_mu + 1/Var_phi2_mu)                    # joint inverse-var
V_wcls   <- V_drwcls   # they coincide here; verified by simulation below

# Oracle d_t^* (since S_t = emptyset, d_t^* is a scalar per t)
d1_star <- -1 / Var_phi1_mu
d2_star <- -1 / Var_phi2_mu

# ----- Storage ---------------------------------------------------------------

beta_unadj  <- numeric(n_reps)
beta_wcls   <- numeric(n_reps)
beta_drwcls <- numeric(n_reps)
beta_opt    <- numeric(n_reps)

# ----- Monte Carlo loop ------------------------------------------------------

cat(sprintf("Running %d Monte Carlo replications, n = %d, alpha = %g, sigma = %g.\n",
            n_reps, n_obs, alpha, sigma))
t0 <- Sys.time()

for (s in seq_len(n_reps)) {
  # Generate one dataset
  A1   <- rbinom(n_obs, 1, 0.5)
  eps1 <- rnorm(n_obs, 0, sigma)
  Y2   <- A1 + eps1

  A2   <- rbinom(n_obs, 1, 0.5)
  eps2 <- rnorm(n_obs, 0, sigma)
  Y3   <- alpha * (Y2 - 0.5) * (A2 - 0.5) + A2 + eps2

  W1 <- 4 * (A1 - 0.5)
  W2 <- 4 * (A2 - 0.5)
  c1 <- A1 - 0.5
  c2 <- A2 - 0.5

  # ---- A. Unadjusted: solve sum_t W_t (Y_{t+1} - (A_t - 0.5) beta) = 0 ----
  num_A <- sum(W1 * Y2 + W2 * Y3)
  den   <- sum(W1 * c1 + W2 * c2)
  beta_unadj[s] <- num_A / den

  # ---- B. WCLS (intercept-only OLS on stacked data) ----
  # WCLS = lm(Y ~ 1 + (A - 0.5)) on n*T rows
  Y_st  <- c(Y2, Y3)
  Ac_st <- c(c1, c2)
  fit   <- lm.fit(x = cbind(1, Ac_st), y = Y_st)
  beta_wcls[s] <- fit$coefficients[2]

  # ---- C. DR-WCLS: our proposed with d_t = 1/4 and mu_t = mu_t^* ----
  # b_t^* = 0.5 in this construction; estimating eq:
  #   sum_t W_t (Y_{t+1} - (A_t-0.5) beta - 0.5) = 0
  num_C <- sum(W1 * (Y2 - 0.5) + W2 * (Y3 - 0.5))
  beta_drwcls[s] <- num_C / den

  # ---- D. Optimal proposed: d_t = d_t^*, mu_t = mu_t^* ----
  # estimating eq:
  #   d1* sum_i W1 (Y2 - (A1-0.5) beta - 0.5) +
  #   d2* sum_i W2 (Y3 - (A2-0.5) beta - 0.5) = 0
  num_D <- d1_star * sum(W1 * (Y2 - 0.5)) + d2_star * sum(W2 * (Y3 - 0.5))
  den_D <- d1_star * sum(W1 * c1)         + d2_star * sum(W2 * c2)
  beta_opt[s] <- num_D / den_D
}

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("Done in %.1f seconds.\n\n", elapsed))

# ----- Diagnostics -----------------------------------------------------------

mc_se_factor <- sqrt(2 / (n_reps - 1))   # SE of n * Var as multiple of n * Var

results <- data.frame(
  estimator = c("A. Unadjusted",
                "B. WCLS (OLS)",
                "C. DR-WCLS",
                "D. Optimal proposed"),
  bias     = c(mean(beta_unadj),
               mean(beta_wcls),
               mean(beta_drwcls),
               mean(beta_opt)) - beta_star,
  nVar_emp = c(n_obs * var(beta_unadj),
               n_obs * var(beta_wcls),
               n_obs * var(beta_drwcls),
               n_obs * var(beta_opt)),
  V_theor  = c(V_unadj, V_wcls, V_drwcls, V_opt)
)
results$mc_se_nVar <- results$nVar_emp * mc_se_factor
results$RE_vs_A    <- results$nVar_emp[1] / results$nVar_emp   # > 1 if more efficient than A

cat("=== Configuration ===\n")
cat(sprintf("  alpha = %g, sigma = %g, n = %d, n_reps = %d\n\n",
            alpha, sigma, n_obs, n_reps))

cat("=== Results ===\n")
print(results, digits = 4, row.names = FALSE)

cat("\n")
cat("Notes:\n")
cat("  bias        = empirical mean(beta_hat) - beta_star (target = 1).\n")
cat("  nVar_emp    = n * empirical Var(beta_hat); converges to V (asymptotic var).\n")
cat("  V_theor     = closed-form asymptotic variance.\n")
cat("  mc_se_nVar  = Monte Carlo standard error of nVar_emp.\n")
cat("  RE_vs_A     = nVar_A / nVar_X, > 1 means more efficient than unadjusted.\n")
cat("\n")
cat("Expected: A and B are within MC error of V_unadj and V_wcls.\n")
cat("          B and C are essentially identical (b_t^* = 0.5 = WCLS-alpha*).\n")
cat("          A is more efficient than B and C: counter-example confirmed for fixed-d.\n")
cat("          A is LESS efficient than D: optimal d compensates for the misadjustment.\n")
cat("\n")

# ----- Optional: alpha sweep, only theoretical curve -------------------------

run_sweep <- TRUE
if (run_sweep) {
  cat("=== Theoretical sweep over alpha ===\n")
  alphas <- c(-5, -4, -3, -2.5, -2, -1.5, -1, 0, 1, 2, 3, 4)
  sweep_df <- data.frame(alpha = alphas)
  for (i in seq_along(alphas)) {
    a <- alphas[i]
    Vp1 <- 4 * sigma^2
    Vp2 <- a^2 * (0.25 + sigma^2) + 4 * sigma^2
    sweep_df$V_unadj[i]    <- (2 + 8*sigma^2 + a^2 * (0.25 + sigma^2) + a) / 4
    sweep_df$V_drwcls[i]   <- (8*sigma^2 + a^2 * (0.25 + sigma^2)) / 4
    sweep_df$V_optimal[i]  <- 1 / (1/Vp1 + 1/Vp2)
    sweep_df$RE_DR_vs_A[i] <- sweep_df$V_unadj[i] / sweep_df$V_drwcls[i]
    sweep_df$RE_Opt_vs_A[i] <- sweep_df$V_unadj[i] / sweep_df$V_optimal[i]
  }
  print(sweep_df, digits = 4, row.names = FALSE)
  cat("\n")
  cat("Reading the sweep:\n")
  cat("  RE_DR_vs_A < 1  means DR-WCLS is LESS efficient than unadjusted.\n")
  cat("  This happens iff alpha < -2 (the counter-example regime).\n")
  cat("  RE_Opt_vs_A > 1 always: the optimal-d proposed estimator is more efficient.\n")
}
