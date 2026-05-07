# simulation_R_code.R
#
# Section 5.4 of the EJS revision: time-varying causal excursion effect.
# Continuous outcome, T = 10, gamma_t(beta) = beta_0 + beta_1 * (t / T).
#
# Four estimators run within each Monte Carlo iteration on the same generated
# dataset (paired Monte Carlo for tighter relative-efficiency comparisons):
#   1. WCLS                          (baseline, no augmentation)
#   2. GAM cross-fitted              (mu_hat fit by mgcv::gam, 10 folds)
#   3. SL.smooth cross-fitted        (mu_hat fit by SL with smooth library)
#   4. Oracle                        (true mu_t, empirical d via fit_d on the truth)
#
# Sweep grid (figures match the layout of Figure 1 in the existing manuscript):
#   left panel:   n in {30, 50, 100} at lambda_1 = 1
#   middle/right: n = 100, lambda_1 in {0, 0.5, ..., 3}
#   crossed with 4 control_patterns (linear, sine, dbeta, step) -> 36 settings
#   1000 MC reps each, packed as nsim_per_seed = 100 x n_seeds = 10 -> 360 tasks
#
# Usage (HPC):  Rscript simulation_R_code.R <itask>

rm(list = ls())

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rootSolve))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(SuperLearner))
suppressPackageStartupMessages(library(earth))

source("function/dgm_cont_time_varying_CEE.R")
source("function/wcls_eif.R")
source("function/eif.R")
source("function/fit_d.R")
source("function/wcls_original.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
    itask <- as.integer(args[1])
} else {
    itask <- 1
    warning("itask set to ", itask, " for debugging")
}

# ---- Simulation design ----------------------------------------------------

# nsim_per_seed can be overridden via the second CLI arg (for smoke testing).
# Production runs use 100, which together with n_seeds = 10 gives 1000 reps.
nsim_per_seed <- if (length(args) >= 2) as.integer(args[2]) else 100
n_seeds <- 10
total_T <- 10

# (sample_size, control_pattern_tuning_param) combos sufficient for Figure 1's
# left/middle/right panels.
ss_lam_grid <- rbind(
    data.frame(sample_size = c(30, 50, 100),
               control_pattern_tuning_param = 1),
    data.frame(sample_size = 100,
               control_pattern_tuning_param = setdiff(seq(0, 3, by = 0.5), 1))
)

# Cross with the 4 mu_t patterns to get 36 settings.
setting_grid <- expand.grid(
    setting_row = 1:nrow(ss_lam_grid),
    control_pattern = c("linear", "sine", "dbeta", "step"),
    stringsAsFactors = FALSE
)
setting_grid <- setting_grid %>%
    mutate(sample_size = ss_lam_grid$sample_size[setting_row],
           control_pattern_tuning_param = ss_lam_grid$control_pattern_tuning_param[setting_row]) %>%
    select(-setting_row) %>%
    mutate(setting_id = row_number())

# Cross settings with seeds to get 360 tasks.
simulation_design <- expand.grid(
    seed = 1:n_seeds,
    setting_id = setting_grid$setting_id,
    stringsAsFactors = FALSE
) %>%
    left_join(setting_grid, by = "setting_id") %>%
    arrange(setting_id, seed)

n_tasks <- nrow(simulation_design)

if (itask < 1 || itask > n_tasks) {
    stop("itask out of range; expected 1..", n_tasks)
}

setting <- simulation_design[itask, ]
sample_size <- setting$sample_size
control_pattern <- setting$control_pattern
control_pattern_tuning_param <- setting$control_pattern_tuning_param
seed <- setting$seed

cat("Conducting simulation for:\n")
print(setting)

set.seed(seed)

# ---- Configurations -------------------------------------------------------

# Same nuisance/error settings as Figure 1 in the existing manuscript.
rand_prob_pattern <- "constant"
rand_prob_tuning_param <- 1
error_var_pattern <- "ind"
error_var_pattern_tuning_param <- 0
error_var_scale_tuning_param <- 1
error_corr_tuning_param <- 0.5

moderator <- "dp_norm"   # the time-moderation column in the data frame
control <- c("dp", "S")  # working-model covariates for fitting mu_t

cf_fold <- 10            # default in wcls_eif

# ---- Helper to run all four estimators on a given dataset -----------------

run_all_methods <- function(dta) {
    out <- list()

    out$wcls <- tryCatch(
        wcls_original(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            numerator_prob = 0.5
        ),
        error = function(e) list(error_message = conditionMessage(e))
    )

    out$gam_cf <- tryCatch(
        wcls_eif(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            ml_method = "gam",
            cross_fit = TRUE,
            cf_fold = cf_fold
        ),
        error = function(e) list(error_message = conditionMessage(e))
    )

    out$sl_smooth_cf <- tryCatch(
        wcls_eif(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            ml_method = "sl.smooth",
            cross_fit = TRUE,
            cf_fold = cf_fold
        ),
        error = function(e) list(error_message = conditionMessage(e))
    )

    # Oracle: use the true mu_t columns expect_Y_A0 / expect_Y_A1 from the dgm,
    # combined with an empirical d_t computed from the truth-evaluated phi^2.
    out$oracle <- tryCatch(
        {
            beta_init <- eif_core(
                dta,
                id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
                moderator = moderator, availability = NULL,
                mu_a0 = "expect_Y_A0", mu_a1 = "expect_Y_A1",
                d_vector = NULL, no_se = TRUE, type = "continuous"
            )$beta_hat
            d_array <- fit_d(
                dta_train = dta, dta_holdout = dta,
                id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
                moderator = moderator, availability = NULL,
                mu_a0 = "expect_Y_A0", mu_a1 = "expect_Y_A1",
                beta = beta_init, d_model_type = "empirical",
                outcome_type = "continuous"
            )
            eif_core(
                dta,
                id = "userid", outcome = "Y", treatment = "A", rand_prob = "prob_A",
                moderator = moderator, availability = NULL,
                mu_a0 = "expect_Y_A0", mu_a1 = "expect_Y_A1",
                d_vector = d_array, no_se = FALSE, type = "continuous"
            )
        },
        error = function(e) list(error_message = conditionMessage(e))
    )

    return(out)
}

# ---- Run the MC iterations for this task ----------------------------------

result_collected <- list()
start_time <- Sys.time()
print_every_n_sims <- 10

for (isim in 1:nsim_per_seed) {
    if (isim %% print_every_n_sims == 0) {
        hours_diff <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
        cat(sprintf("Starting isim: %d/%d; Hours lapsed: %s\n",
                    isim, nsim_per_seed, hours_diff))
    }

    dta <- dgm_cont_time_varying_CEE(
        sample_size = sample_size,
        total_T = total_T,
        rand_prob_pattern = rand_prob_pattern,
        rand_prob_tuning_param = rand_prob_tuning_param,
        control_pattern = control_pattern,
        control_pattern_tuning_param = control_pattern_tuning_param,
        error_var_pattern = error_var_pattern,
        error_var_pattern_tuning_param = error_var_pattern_tuning_param,
        error_var_scale_tuning_param = error_var_scale_tuning_param,
        error_corr_tuning_param = error_corr_tuning_param
    )

    result_collected[[isim]] <- run_all_methods(dta)
}

cat("Last iteration result preview:\n")
str(result_collected[[length(result_collected)]], max.level = 2)

dir.create("result_tmp", showWarnings = FALSE)
saveRDS(result_collected, file = paste0("result_tmp/", itask, ".RDS"))

# Result aggregation is handled by the standalone collect_results.R.
