rm(list = ls())

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
jobid <- as.integer(args[1])

nsim <- 1000

source("dgm_cont_time_s_tvvar.R")

set.seed(123)

# varying control_pattern
SD1 <- expand.grid(sample_size = c(30, 50, 100),
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
#                    stringsAsFactors = FALSE)
# 
# SD <- SD1

# varying error_var_pattern
SD2 <- expand.grid(sample_size = c(30, 50, 100),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = c(1,2),
                   error_var_pattern = "timevar-linear",
                   error_var_pattern_tuning_param = seq(from = 0, to = 3, by = 0.5),
                   error_var_scale_tuning_param = 1,
                   error_corr_tuning_param = 0.5,
                   stringsAsFactors = FALSE)

# varying error_var_scale
SD3 <- expand.grid(sample_size = c(30, 50, 100),
                   total_T = 10,
                   S_pattern = c("constant"),
                   rand_prob_pattern = c("constant"),
                   rand_prob_tuning_param = 1,
                   control_pattern = c("linear", "sine", "dbeta", "step"),
                   control_pattern_tuning_param = c(1,2),
                   error_var_pattern = "ind",
                   error_var_pattern_tuning_param = 0,
                   error_var_scale_tuning_param = seq(from = 1, to = 5, by = 1),
                   error_corr_tuning_param = 0.5,
                   stringsAsFactors = FALSE)

SD <- rbind(SD1, SD2, SD3)

SD <- SD %>% distinct()

# for (jobid in 1:nrow(SD)) {

print(SD[jobid, ])

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
                   "-nsim=", nsim,
                   ".RDS")

datasets <- list()
for (isim in 1:nsim) {
    dta <- dgm_cont_time_s_tvvar(
        sample_size,
        total_T,
        S_pattern,
        rand_prob_pattern,
        rand_prob_tuning_param,
        control_pattern,
        control_pattern_tuning_param,
        error_var_pattern,
        error_var_pattern_tuning_param,
        error_var_scale_tuning_param,
        error_corr_tuning_param
    )
    datasets <- c(datasets, list(dta))
}

dir.create("dgm_cont_time_s_tvvar", showWarnings = FALSE)
saveRDS(datasets, file = paste0("dgm_cont_time_s_tvvar/", filename))    
# }
