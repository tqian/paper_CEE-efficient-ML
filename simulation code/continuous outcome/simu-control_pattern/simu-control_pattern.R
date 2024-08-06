# continuous outcome
rm(list = ls())

nsim_overall <- 1000

args <- commandArgs(trailingOnly = TRUE)
jobid <- as.integer(args[1])

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rootSolve))

source("../functions/wcls_eif.R")
source("../functions/fit_d.R")
source("../functions/eif.R")
source("../functions/wcls_original.R")

# 1. Simulation design setup ----------------------------------------------

## fast: gam, lm, oracle, wcls
SD1 <- expand.grid(sample_size = c(50),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                   error_var_pattern = "ind",
                   error_var_pattern_tuning_param = 0,
                   error_var_scale_tuning_param = 1,
                   error_corr_tuning_param = 0.5,
                   moderator = c("NULL"),
                   control = c("dp", "dp S"),
                   estimator = "eif",
                   ml_method = c("gam", "lm"),
                   cross_fit = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)

# SD1 <- expand.grid(sample_size = c(50),
#                    total_T = 10,
#                    S_pattern = c("constant"),
#                    rand_prob_pattern = c("constant"),
#                    rand_prob_tuning_param = 1,
#                    control_pattern = c("linear"),
#                    control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
#                    error_var_pattern = "ind",
#                    error_var_pattern_tuning_param = 0,
#                    error_var_scale_tuning_param = 1,
#                    error_corr_tuning_param = 0.5,
#                    moderator = c("NULL"),
#                    control = c("dp", "dp S"),
#                    estimator = "eif",
#                    ml_method = c("gam", "lm"),
#                    cross_fit = c(FALSE),
#                    stringsAsFactors = FALSE)
# 
# SD <- SD1
# SD <- SD %>% mutate(subtaskid = 1, subtasktotal = 1, taskspeed = "fast")

SD2 <- expand.grid(sample_size = c(50),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                   error_var_pattern = "ind",
                   error_var_pattern_tuning_param = 0,
                   error_var_scale_tuning_param = 1,
                   error_corr_tuning_param = 0.5,
                   moderator = c("NULL"),
                   control = c("dp", "dp S"),
                   estimator = c("eif_oracle", "wcls"),
                   ml_method = "",
                   cross_fit = FALSE,
                   stringsAsFactors = FALSE)

SD_fast <- rbind(SD1, SD2)
SD_fast <- SD_fast %>% mutate(subtaskid = 1, subtasktotal = 1, taskspeed = "fast")

## slower: ranger, sl.smooth
SD_slower1 <- expand.grid(sample_size = c(50),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "ranger",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:10,
                          subtasktotal = 10,
                          taskspeed = "slower",
                          stringsAsFactors = FALSE)

SD_slower2 <- expand.grid(sample_size = c(50),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "sl.smooth",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:10,
                          subtasktotal = 10,
                          taskspeed = "slower",
                          stringsAsFactors = FALSE)

SD_slower <- rbind(SD_slower1, SD_slower2)

## slowest: sl.all
SD_slowest <- expand.grid(sample_size = c(50),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "sl.all",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:100,
                          subtasktotal = 100,
                          taskspeed = "slowest",
                          stringsAsFactors = FALSE)

SD <- rbind(SD_fast, SD_slower, SD_slowest)


# 1.1. Additional simulation design setup -------------------------------------


## fast: gam, lm, oracle, wcls
SD1 <- expand.grid(sample_size = c(30, 100),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = 1,
                   error_var_pattern = "ind",
                   error_var_pattern_tuning_param = 0,
                   error_var_scale_tuning_param = 1,
                   error_corr_tuning_param = 0.5,
                   moderator = c("NULL"),
                   control = c("dp", "dp S"),
                   estimator = "eif",
                   ml_method = c("gam", "lm"),
                   cross_fit = c(TRUE, FALSE),
                   stringsAsFactors = FALSE)

SD2 <- expand.grid(sample_size = c(30, 100),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = 1,
                   error_var_pattern = "ind",
                   error_var_pattern_tuning_param = 0,
                   error_var_scale_tuning_param = 1,
                   error_corr_tuning_param = 0.5,
                   moderator = c("NULL"),
                   control = c("dp", "dp S"),
                   estimator = c("eif_oracle", "wcls"),
                   ml_method = "",
                   cross_fit = FALSE,
                   stringsAsFactors = FALSE)

SD_fast <- rbind(SD1, SD2)
SD_fast <- SD_fast %>% mutate(subtaskid = 1, subtasktotal = 1, taskspeed = "fast")

## slower: ranger, sl.smooth
SD_slower1 <- expand.grid(sample_size = c(30, 100),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = 1,
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "ranger",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:10,
                          subtasktotal = 10,
                          taskspeed = "slower",
                          stringsAsFactors = FALSE)

SD_slower2 <- expand.grid(sample_size = c(30, 100),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = 1,
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "sl.smooth",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:10,
                          subtasktotal = 10,
                          taskspeed = "slower",
                          stringsAsFactors = FALSE)

SD_slower <- rbind(SD_slower1, SD_slower2)

## slowest: sl.all
SD_slowest <- expand.grid(sample_size = c(30, 100),
                          total_T = 10,
                          S_pattern = c("constant"),
                          rand_prob_pattern = c("constant"),
                          rand_prob_tuning_param = 1,
                          control_pattern = c("linear", "sine", "dbeta", "step"),
                          control_pattern_tuning_param = 1,
                          error_var_pattern = "ind",
                          error_var_pattern_tuning_param = 0,
                          error_var_scale_tuning_param = 1,
                          error_corr_tuning_param = 0.5,
                          moderator = c("NULL"),
                          control = c("dp", "dp S"),
                          estimator = c("eif"),
                          ml_method = "sl.all",
                          cross_fit = c(TRUE, FALSE),
                          subtaskid = 1:100,
                          subtasktotal = 100,
                          taskspeed = "slowest",
                          stringsAsFactors = FALSE)

SD2 <- rbind(SD_fast, SD_slower, SD_slowest)

SD <- rbind(SD, SD2)

# 1.2. Configure the current job ------------------------------------------

for (jobid in 1:nrow(SD)) {
  
print(SD[jobid, ])

if (jobid == 1) {
    write_csv(SD, file = "simulation_design.csv")
}


## simulation setting specification

S_pattern <- SD$S_pattern[jobid]
rand_prob_pattern <- SD$rand_prob_pattern[jobid]
rand_prob_tuning_param <- SD$rand_prob_tuning_param[jobid]
control_pattern <- SD$control_pattern[jobid]
control_pattern_tuning_param <- SD$control_pattern_tuning_param[jobid]
error_var_pattern <- SD$error_var_pattern[jobid]
error_var_pattern_tuning_param <- SD$error_var_pattern_tuning_param[jobid]
error_var_scale_tuning_param <- SD$error_var_scale_tuning_param[jobid]
error_corr_tuning_param <- SD$error_corr_tuning_param[jobid]
sample_size <- SD$sample_size[jobid]
total_T <- SD$total_T[jobid]

moderator <- SD$moderator[jobid]
control <- SD$control[jobid]
estimator <- SD$estimator[jobid]
ml_method <- SD$ml_method[jobid]
cross_fit <- SD$cross_fit[jobid]

d_true <- -1 * (4 * (0.2^2 / 3 +
                         seq(from = error_var_scale_tuning_param,
                             by = error_var_pattern_tuning_param,
                             length = total_T)))^(-1)

if (ml_method %in% c("sl.smooth", "sl.all")) {
    suppressPackageStartupMessages(library(SuperLearner))
    suppressPackageStartupMessages(library(earth))
    suppressPackageStartupMessages(library(ranger))
    suppressPackageStartupMessages(library(nnet))
    suppressPackageStartupMessages(library(xgboost))
} else if (ml_method == "ranger") {
    suppressPackageStartupMessages(library(ranger))
} else if (ml_method == "gam") {
    suppressPackageStartupMessages(library(mgcv))
} else if (ml_method == "rf") {
    suppressPackageStartupMessages(library(randomForest))
}

if (moderator == "NULL") {
    moderator <- NULL
}

control <- unlist(strsplit(control, " "))


## parallel job specification

subtaskid <- SD$subtaskid[jobid]
subtasktotal <- SD$subtasktotal[jobid]
taskspeed <- SD$taskspeed[jobid]

print_every_n_sims <- case_when(taskspeed == "fast" ~ 100,
                                taskspeed == "slower" ~ 10,
                                taskspeed == "slowest" ~ 1)

nsim_this_job <- nsim_overall / subtasktotal
isim_start <- nsim_this_job * (subtaskid - 1) + 1
isim_end <- nsim_this_job * subtaskid

set.seed(subtaskid)


# 2. Run simulations ------------------------------------------------------


filename <- paste0("dgm_cont_time_s_tvvar",
                   "-S=", S_pattern,
                   "-rand=", rand_prob_pattern,
                   "-rand_param=", rand_prob_tuning_param,
                   "-ctrl=", control_pattern,
                   "-ctrl_param=", control_pattern_tuning_param,
                   "-errvar=", error_var_pattern,
                   "-errvar_param=", error_var_pattern_tuning_param,
                   "-errvarscale_param=", error_var_scale_tuning_param,
                   "errcorr_param=", error_corr_tuning_param,
                   "-ss=", sample_size,
                   "-T=", total_T,
                   "-nsim=", 1000,
                   ".RDS")

#filename <- "dgm_cont_time_s_tvvar-S=constant-rand=constant-rand_param=1-ctrl=sine-ctrl_param=1-errvar=timevar-linear-errvar_param=3-errvarscale_param=1errcorr_param=0.5-ss=50-T=10-nsim=1000.RDS"

datasets <- readRDS(paste0("../datasets/dgm_cont_time_s_tvvar/", filename))

start_time <- Sys.time()

fit_list <- c()
for (isim in isim_start:isim_end) {
    if (isim %% print_every_n_sims == 0) {
        current_time <- Sys.time()
        hours_diff <- round(difftime(current_time, start_time, units = "hours"), 2)
        cat(paste0("Starting isim: ", isim, "/", nsim_this_job, "; Hours lapsed: ", hours_diff, "\n"))
    }
    dta <- datasets[[isim]]

    if (estimator == "eif") {

        fit <- wcls_eif(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            ml_method = ml_method,
            cross_fit = cross_fit
        )
    } else if (estimator == "eif_oracle") {
        fit <- wcls_eif_oracle(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            mu_a0 = "expect_Y_A0",
            mu_a1 = "expect_Y_A1",
            d_true = d_true
        )
    } else if (estimator == "wcls") {
        fit <- wcls_original(
            dta = dta,
            id = "userid",
            outcome = "Y",
            treatment = "A",
            rand_prob = "prob_A",
            moderator = moderator,
            control = control,
            numerator_prob = 0.5
        )
    }

    fit_list <- c(fit_list, list(fit))
}

dir.create("result", showWarnings = FALSE)

print("The last model fit is:")
print(fit)

saveRDS(fit_list, file = paste0("result/", jobid, ".RDS"))
}


# Collect all results and create a big data frame -------------------------

if (jobid == nrow(SD)) {
    
    print(paste0("Collecting ", nrow(SD), " result files..."))
    
    all_result_files <- paste0("result/", 1:nrow(SD), ".RDS")
    
    while (!all(file.exists(all_result_files))) {
        Sys.sleep(5)
    }
    
    result_all <- data.frame()
    
    beta0_true <- 0.5
    beta1_true <- 0.2
    
    for (jobid in 1:nrow(SD)) {
        if (jobid %% 10 == 0) {
            cat(jobid, "")
        }
        
        fit_list <- readRDS(paste0("result/", jobid, ".RDS"))
        
        if (SD$moderator[jobid] == "NULL") {
            # only beta0
            
            beta0 <- sapply(fit_list, function(l) l$beta_hat)
            # use adjusted se instead
            beta0_se <- sapply(fit_list, function(l) l$beta_se)
            beta0_lci <- beta0 - 1.96 * beta0_se
            beta0_uci <- beta0 + 1.96 * beta0_se
            
            beta1 <- NA
            beta1_se <- NA
            beta1_lci <- NA
            beta1_uci <- NA
            
        } else if (SD$moderator[jobid] == "S") {
            # both beta0 and beta1
            stop("Haven't revised the following code for eif simulation yet")
            
            est <- t(sapply(fit_list, function(l) l$beta_hat))
            se <- t(sapply(fit_list, function(l) l$beta_se_adjusted))
            lci <- t(sapply(fit_list, function(l) l$conf_int_adjusted[, 1]))
            uci <- t(sapply(fit_list, function(l) l$conf_int_adjusted[, 2]))
            
            beta0 <- est[, 1]
            beta0_se <- se[, 1]
            beta0_lci <- lci[, 1]
            beta0_uci <- uci[, 1]
            
            beta1 <- est[, 2]
            beta1_se <- se[, 2]
            beta1_lci <- lci[, 2]
            beta1_uci <- uci[, 2]
        }
        
        result_df <- data.frame(
            outcome = "continuous",
            S_pattern = SD$S_pattern[jobid],
            rand_prob_pattern = SD$rand_prob_pattern[jobid],
            rand_prob_tuning_param = SD$rand_prob_tuning_param[jobid],
            control_pattern = SD$control_pattern[jobid],
            control_pattern_tuning_param = SD$control_pattern_tuning_param[jobid],
            error_var_pattern = SD$error_var_pattern[jobid],
            error_var_pattern_tuning_param = SD$error_var_pattern_tuning_param[jobid],
            error_var_scale_tuning_param = SD$error_var_scale_tuning_param[jobid],
            error_corr_tuning_param = SD$error_corr_tuning_param[jobid],
            sample_size = SD$sample_size[jobid],
            total_T = SD$total_T[jobid],
            moderator = SD$moderator[jobid],
            control = SD$control[jobid],
            estimator = SD$estimator[jobid],
            ml_method = SD$ml_method[jobid],
            cross_fit = SD$cross_fit[jobid],
            beta0 = beta0,
            beta0_se = beta0_se,
            beta0_lci = beta0_lci,
            beta0_uci = beta0_uci,
            beta1 = beta1,
            beta1_se = beta1_se,
            beta1_lci = beta1_lci,
            beta1_uci = beta1_uci
        )
        
        result_all <- rbind(result_all, result_df)
    }
    
    saveRDS(result_all, "result_all.RDS")
}
