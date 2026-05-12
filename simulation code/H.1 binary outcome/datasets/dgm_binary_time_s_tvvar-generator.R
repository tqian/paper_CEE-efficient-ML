rm(list = ls())

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
jobid <- as.integer(args[1])

nsim <- 1000

source("dgm_binary_time_s_tvvar.R")

set.seed(123)

# varying control_pattern -- can't separate error variance and mean
SD1 <- expand.grid(sample_size = c(30, 50, 100),
                  total_T = 10,
                  S_pattern = c("constant"),
                  rand_prob_pattern = c("constant"),
                  rand_prob_tuning_param = 1,
                  control_pattern = c("linear", "sine", "dbeta", "step"),
                  control_tuning_param = seq(from = 0, to = 1, by = 0.2),
                  ar1_tuning_param = 0.1,
                  stringsAsFactors = FALSE)

SD <- SD1

SD <- SD %>% distinct()

for (jobid in 1:nrow(SD)) {
print(SD[jobid, ])

S_pattern <- SD$S_pattern[jobid]
rand_prob_pattern <- SD$rand_prob_pattern[jobid]
rand_prob_tuning_param <- SD$rand_prob_tuning_param[jobid]
control_pattern <- SD$control_pattern[jobid]
control_pattern_tuning_param <- SD$control_tuning_param[jobid]
ar1_tuning_param <- SD$ar1_tuning_param[jobid]
sample_size <- SD$sample_size[jobid]
total_T <- SD$total_T[jobid]

filename <- paste0("dgm_binary_time_s_tvvar",
                   "-S=", S_pattern,
                   "-rand=", rand_prob_pattern,
                   "-rand_param=", rand_prob_tuning_param,
                   "-ctrl=", control_pattern,
                   "-ctrl_param=", control_pattern_tuning_param,
                   "-ar1_param=", ar1_tuning_param,
                   "-ss=", sample_size,
                   "-T=", total_T,
                   "-nsim=", nsim,
                   ".RDS")

datasets <- list()
for (isim in 1:nsim) {
    dta <- dgm_binary_time_s_tvvar(
      sample_size = sample_size, 
      total_T = total_T,
      S_pattern = S_pattern,
      rand_prob_pattern = rand_prob_pattern,
      rand_prob_tuning_param = rand_prob_tuning_param,
      control_pattern = control_pattern,
      control_tuning_param = control_pattern_tuning_param, 
      ar1_tuning_param = ar1_tuning_param)

    datasets <- c(datasets, list(dta))
}


dir.create("dgm_binary_time_s_tvvar", showWarnings = FALSE)
saveRDS(datasets, file = paste0("dgm_binary_time_s_tvvar/", filename))   
}
    
