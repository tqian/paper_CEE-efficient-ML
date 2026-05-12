# =============================================================================
# Sharper counter-example: even the OPTIMAL proposed estimator (d_t^*, mu_t^*)
# loses to the unadjusted estimator.
#
# DGP:
#   T = 2, p_t = 0.5, I_t = 1.
#   Y_2 = A_1 + eps_1
#   Y_3 = theta * (A_1 - 0.5) * (A_2 - 0.5) + A_2 + eps_2,   theta = -3.
#   Working CEE gamma_t(emptyset; beta) = beta. beta_star = 1.
#
# Estimators (all called via the project's own functions):
#   A. Unadjusted        : eif_core(...) with mu_hat_a0 = mu_hat_a1 = 0,
#                                            d_vector = NULL  (defaults to 1).
#   B. WCLS              : wcls_original(...) with control = NULL (intercept-only).
#   C. DR-WCLS           : eif_core(...) with mu_hat_a* = mu_t^*,
#                                            d_vector = NULL  (defaults to 1).
#   D. Optimal proposed  : eif_core(...) with mu_hat_a* = mu_t^*,
#                                            d_vector = oracle d_t^*.
#
# Theory (with theta = -3, sigma^2 = 1; beta_star = 1; A_bread = -2; V = B / 4):
#   V_unadj   = 2.3125
#   V_DRWCLS  = 2.5625   (= V_WCLS in this construction since b_t^* = WCLS-alpha*)
#   V_optimal = 2.4390
# All three covariate-adjusted estimators are LESS efficient than unadjusted.
# =============================================================================

library(rootSolve)

# ----- Source the project's own estimator functions --------------------------
funcs_dir <- "../F.1 wallclock benchmark/functions"
source(file.path(funcs_dir, "eif.R"))
source(file.path(funcs_dir, "wcls_original.R"))

# ----- Configuration ---------------------------------------------------------

theta   <- -3
sigma   <- 1
n_obs   <- 2000
n_reps  <- 5000
T_obs   <- 2
beta_star <- 1

# ----- Theoretical reference values ------------------------------------------

Var_phi1_mu <- 4 * sigma^2                        # = 4
Var_phi2_mu <- 0.25 * theta^2 + 4 * sigma^2       # = 6.25
Cov_mu_star <- 0                                  # (WA-2) at mu* holds

V_unadj   <- ((1 + 4*sigma^2) + (0.25*theta^2 + 1 + 4*sigma^2) + 2 * 0.5 * theta) / 4
V_drwcls  <- (Var_phi1_mu + Var_phi2_mu + 2 * Cov_mu_star) / 4
V_opt     <- (Var_phi1_mu * Var_phi2_mu - Cov_mu_star^2) /
             (Var_phi1_mu + Var_phi2_mu - 2 * Cov_mu_star)
V_wcls    <- V_drwcls

# Oracle d_t^* (S_t = empty, so it's a single number per t)
d_oracle_per_t <- c(-1 / Var_phi1_mu, -1 / Var_phi2_mu)

# ----- Data generation: long format expected by eif_core / wcls_original -----
# Rows are (id, t) pairs, sorted by id then t.

simulate_dataset <- function(n, theta, sigma) {
  A1   <- rbinom(n, 1, 0.5)
  eps1 <- rnorm(n,  0, sigma)
  Y2   <- A1 + eps1

  A2   <- rbinom(n, 1, 0.5)
  eps2 <- rnorm(n,  0, sigma)
  Y3   <- theta * (A1 - 0.5) * (A2 - 0.5) + A2 + eps2

  # True conditional means mu_t^*(H_t, A_t = a) = E[Y_{t+1} | H_t, A_t = a].
  # t = 1: H_1 = empty, mu_1^*(., 0) = 0, mu_1^*(., 1) = 1.
  mu1_a0 <- rep(0, n)
  mu1_a1 <- rep(1, n)
  # t = 2: mu_2^*(H_2, a) = theta (A_1 - 0.5)(a - 0.5) + a.
  mu2_a0 <- -0.5 * theta * (A1 - 0.5)
  mu2_a1 <-  0.5 * theta * (A1 - 0.5) + 1

  data.frame(
    userid    = rep(1:n, each = T_obs),
    t         = rep(1:T_obs, times = n),
    A         = as.vector(rbind(A1, A2)),
    Y         = as.vector(rbind(Y2, Y3)),
    prob_A    = 0.5,
    mu_hat_a0 = as.vector(rbind(mu1_a0, mu2_a0)),
    mu_hat_a1 = as.vector(rbind(mu1_a1, mu2_a1))
  )
}

# ----- Wrappers around the project functions ---------------------------------

fit_unadjusted <- function(dta) {
  d2 <- dta
  d2$mu_hat_a0 <- 0
  d2$mu_hat_a1 <- 0
  fit <- eif_core(d2, id = "userid", outcome = "Y", treatment = "A",
                  rand_prob = "prob_A", moderator = NULL, availability = NULL,
                  mu_a0 = "mu_hat_a0", mu_a1 = "mu_hat_a1",
                  d_vector = NULL, no_se = TRUE, type = "continuous")
  as.numeric(fit$beta_hat)
}

fit_drwcls <- function(dta) {
  fit <- eif_core(dta, id = "userid", outcome = "Y", treatment = "A",
                  rand_prob = "prob_A", moderator = NULL, availability = NULL,
                  mu_a0 = "mu_hat_a0", mu_a1 = "mu_hat_a1",
                  d_vector = NULL, no_se = TRUE, type = "continuous")
  as.numeric(fit$beta_hat)
}

fit_optimal <- function(dta, d_oracle_per_t) {
  d_array <- rep(d_oracle_per_t, times = nrow(dta) / length(d_oracle_per_t))
  fit <- eif_core(dta, id = "userid", outcome = "Y", treatment = "A",
                  rand_prob = "prob_A", moderator = NULL, availability = NULL,
                  mu_a0 = "mu_hat_a0", mu_a1 = "mu_hat_a1",
                  d_vector = d_array, no_se = TRUE, type = "continuous")
  as.numeric(fit$beta_hat)
}

fit_wcls <- function(dta) {
  fit <- wcls_original(dta, id = "userid", outcome = "Y", treatment = "A",
                       rand_prob = "prob_A", moderator = NULL, control = NULL,
                       availability = NULL, numerator_prob = 0.5,
                       no_se = TRUE)
  as.numeric(fit$beta_hat)
}

# ----- Monte Carlo loop ------------------------------------------------------

set.seed(20260501)

beta_unadj  <- numeric(n_reps)
beta_wcls   <- numeric(n_reps)
beta_drwcls <- numeric(n_reps)
beta_opt    <- numeric(n_reps)

cat(sprintf("Running %d Monte Carlo replications, n = %d, theta = %g.\n",
            n_reps, n_obs, theta))
t0 <- Sys.time()

for (s in seq_len(n_reps)) {
  dta <- simulate_dataset(n_obs, theta, sigma)
  beta_unadj[s]  <- fit_unadjusted(dta)
  beta_wcls[s]   <- fit_wcls(dta)
  beta_drwcls[s] <- fit_drwcls(dta)
  beta_opt[s]    <- fit_optimal(dta, d_oracle_per_t)
}

elapsed <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("Done in %.1f seconds.\n\n", elapsed))

# ----- Report ----------------------------------------------------------------

mc_se_factor <- sqrt(2 / (n_reps - 1))

results <- data.frame(
  estimator = c("A. Unadjusted (eif_core, mu = 0)",
                "B. WCLS (wcls_original)",
                "C. DR-WCLS (eif_core, mu = mu*, d = 1)",
                "D. Optimal (eif_core, mu = mu*, d = d*)"),
  bias     = c(mean(beta_unadj), mean(beta_wcls),
               mean(beta_drwcls), mean(beta_opt)) - beta_star,
  nVar_emp = c(n_obs * var(beta_unadj),
               n_obs * var(beta_wcls),
               n_obs * var(beta_drwcls),
               n_obs * var(beta_opt)),
  V_theor  = c(V_unadj, V_wcls, V_drwcls, V_opt)
)
results$mc_se_nVar <- results$nVar_emp * mc_se_factor
results$RE_vs_A    <- results$nVar_emp[1] / results$nVar_emp

cat("=== Configuration ===\n")
cat(sprintf("  theta = %g, sigma = %g, n = %d, n_reps = %d\n\n",
            theta, sigma, n_obs, n_reps))

cat("=== Results ===\n")
print(results, digits = 4, row.names = FALSE)

cat("\n")
cat("Counter-example confirmed if RE_vs_A < 1 for D (optimal proposed).\n")
