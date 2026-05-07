# Wall-clock benchmark of nuisance-fit methods on Drink Less data.
#
# Added 2026.05 in response to Reviewer 1 comment 3 (computational scalability).
# Times Algorithm 1 (no cross-fit) and Algorithm 2 (cross-fit, K = 5) for the
# continuous proximal outcome on the Drink Less MRT dataset (n = 349, T = 30).
# Single thread, single run, on the local machine.

rm(list = ls())

library(tidyverse)
library(rootSolve)
library(SuperLearner)
library(mgcv)
library(ranger)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("functions/wcls_original.R")
source("functions/wcls_eif.R")
source("functions/eif.R")
source("functions/fit_d.R")


# 1. Data preparation -----------------------------------------------------

dta <- readRDS("data/DrinkLess_cont_bin.RDS")
dta$ID <- as.factor(dta$ID)

dta <- dta %>%
    group_by(ID) %>%
    mutate(treatment_lag1 = dplyr::lag(treatment, default = 0),
           binary_8to759nextday_lag1 = dplyr::lag(binary_8to759nextday, default = 0),
           seconds_8to759nextday_lag1 = dplyr::lag(seconds_8to759nextday, default = 0),
           seconds_8to9_lag1 = dplyr::lag(seconds_8to9, default = 0),
           binary_8to9_lag1 = dplyr::lag(binary_8to9, default = 0),
           seconds_8to759nextday_lag2 = dplyr::lag(seconds_8to759nextday, default = 0, n = 2)) %>%
    ungroup()
dta <- as.data.frame(dta)
dta$userid <- as.numeric(factor(dta$ID))
dta <- dta[order(dta$userid, dta$decision_index, decreasing = FALSE), ]

n_id <- length(unique(dta$userid))
n_T  <- nrow(dta) / n_id
cat(sprintf("Drink Less data: n = %d participants, T = %d decision points, total rows = %d\n",
            n_id, n_T, nrow(dta)))


# 2. Fixed analysis spec --------------------------------------------------

id <- "userid"
outcome <- "seconds_8to9"
treatment <- "treatment"
rand_prob <- "rand_prob"
moderator <- NULL  # marginal CEE
control <- c("age", "AUDIT_score", "decision_index",
             "seconds_8to759nextday_lag1", "treatment_lag1",
             "seconds_8to759nextday_lag2")
gam_control_spline_var <- c("age", "AUDIT_score", "decision_index")
d_model <- "earth"


# 3. Benchmark wrapper ----------------------------------------------------

time_one <- function(label, expr_quoted) {
    cat(sprintf("\n--- %s ---\n", label))
    set.seed(1)
    t0 <- Sys.time()
    res <- tryCatch(eval(expr_quoted),
                    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
    t1 <- Sys.time()
    elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
    cat(sprintf("Wall clock: %.2f s\n", elapsed))
    list(label = label, elapsed_sec = elapsed, ok = !is.null(res))
}


# 4. Run all methods ------------------------------------------------------

results <- list()

results$wcls <- time_one("WCLS (baseline)",
    quote(wcls_original(dta, id = id, outcome = outcome, treatment = treatment,
                        rand_prob = rand_prob, moderator = moderator, control = control,
                        numerator_prob = 0.6)))

results$gam <- time_one("GAM (no cross-fit)",
    quote(wcls_eif(dta = dta, id = id, outcome = outcome, treatment = treatment,
                   rand_prob = rand_prob, moderator = moderator, control = control,
                   ml_method = "gam", cross_fit = FALSE, d_model_type = d_model,
                   gam_control_spline_var = gam_control_spline_var)))

results$gam_cf <- time_one("GAM (cross-fit, K=5)",
    quote(wcls_eif(dta = dta, id = id, outcome = outcome, treatment = treatment,
                   rand_prob = rand_prob, moderator = moderator, control = control,
                   ml_method = "gam", cross_fit = TRUE, cf_fold = 5,
                   d_model_type = d_model,
                   gam_control_spline_var = gam_control_spline_var)))

results$rf_cf <- time_one("RF / ranger (cross-fit, K=5)",
    quote(wcls_eif(dta = dta, id = id, outcome = outcome, treatment = treatment,
                   rand_prob = rand_prob, moderator = moderator, control = control,
                   ml_method = "ranger", cross_fit = TRUE, cf_fold = 5,
                   d_model_type = d_model,
                   gam_control_spline_var = gam_control_spline_var)))

results$sl_smooth_cf <- time_one("SL.smooth (cross-fit, K=5)  [GLM+GAM+earth, no XGBoost/NN/RF]",
    quote(wcls_eif(dta = dta, id = id, outcome = outcome, treatment = treatment,
                   rand_prob = rand_prob, moderator = moderator, control = control,
                   ml_method = "sl.smooth", cross_fit = TRUE, cf_fold = 5,
                   d_model_type = d_model,
                   gam_control_spline_var = gam_control_spline_var)))

results$sl_all_cf <- time_one("SL.all (cross-fit, K=5)  [includes XGBoost + NN]",
    quote(wcls_eif(dta = dta, id = id, outcome = outcome, treatment = treatment,
                   rand_prob = rand_prob, moderator = moderator, control = control,
                   ml_method = "sl.all", cross_fit = TRUE, cf_fold = 5,
                   d_model_type = d_model,
                   gam_control_spline_var = gam_control_spline_var)))


# 5. Save and print summary -----------------------------------------------

dir.create("result", showWarnings = FALSE)

bench_df <- do.call(rbind, lapply(results, function(r)
    data.frame(method = r$label, elapsed_sec = r$elapsed_sec, ok = r$ok)))
rownames(bench_df) <- NULL

saveRDS(bench_df, "result/wallclock_continuous.RDS")
write.csv(bench_df, "result/wallclock_continuous.csv", row.names = FALSE)

cat("\n\n========== Wall-clock summary (continuous outcome) ==========\n")
print(bench_df, row.names = FALSE)
cat(sprintf("\nR version: %s, machine: %s\n", R.version.string, Sys.info()["nodename"]))
cat(sprintf("Time stamp: %s\n", format(Sys.time())))
