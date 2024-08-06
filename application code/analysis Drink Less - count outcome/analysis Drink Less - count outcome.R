rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(rootSolve)
library(latex2exp)
library(SuperLearner)
library(mgcv)
library(ranger)
library(patchwork) # to arrange plots

library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

source("functions/wcls_original_binary.R")
source("functions/eif.R")
source("functions/emee_eif.R")
source("functions/fit_d.R")
# source("functions/wcls_original_binary.R")
# source("functions/wcls_ml_binary.R")


# 1. preparation ----------------------------------------------------------

dta <- readRDS("data/DrinkLess_cont_bin_count.RDS")
dta$ID <- as.factor(dta$ID)

dta <- dta %>% 
  group_by(ID) %>% 
  mutate(treatment_lag1 = dplyr::lag(treatment, default = 0),
         binary_8to759nextday_lag1 = dplyr::lag(binary_8to759nextday, default = 0),
         seconds_8to759nextday_lag1 = dplyr::lag(seconds_8to759nextday, default = 0),
         count_8to759nextday_lag1 = dplyr::lag(count_8to759nextday, default = 0),
         seconds_8to9_lag1 = dplyr::lag(seconds_8to9, default = 0),
         count_8to9_lag1 = dplyr::lag(count_8to9, default = 0),
         binary_8to9_lag1 = dplyr::lag(binary_8to9, default = 0),
         seconds_8to759nextday_lag2 = dplyr::lag(seconds_8to759nextday, default = 0, n = 2),
         count_8to759nextday_lag2 = dplyr::lag(count_8to759nextday, default = 0, n = 2)) %>% 
  ungroup()
dta <- as.data.frame(dta)

dta$userid <- as.numeric(factor(dta$ID))


# 2. Continuous outcome analysis ------------------------------------------

id <- "userid"
treatment <- "treatment"
rand_prob <- "rand_prob"
control <- c("age", "AUDIT_score", "decision_index", "count_8to759nextday_lag1", "seconds_8to759nextday_lag1",
             "treatment_lag1", "count_8to759nextday_lag2", "seconds_8to759nextday_lag2")
gam_control_spline_var <- c("age", "AUDIT_score", "decision_index")

dta <- dta[order(dta$userid, dta$decision_index, decreasing = F), ]

design_count <- expand.grid(moderator = c("NULL", "decision_index"),
                           outcome = c("count_8to9"),
                           stringsAsFactors = FALSE)

set.seed(1)

for (idesign in 1:nrow(design_count)) {
  
  print(paste0("idesign = ", idesign))
  
  outcome <- design_count$outcome[idesign]
  moderator <- design_count$moderator[idesign]
  model_for_var <- design_count$model_for_var[idesign]
  if (moderator == "NULL") {
    d_model <- "earth"
  } else {
    d_model <- "earth"
  }
  
  if (moderator == "NULL") {
    moderator <- NULL
  }
  
  fit_emee <- wcls_original_binary(
    dta,
    id = id,
    outcome = outcome,
    treatment = treatment,
    rand_prob = rand_prob,
    moderator = moderator,
    control = control,
    numerator_prob = 0.6,
  )
  
  fit_gam <- emee_eif(
    dta = dta,
    id = id,
    outcome = outcome,
    treatment = treatment,
    rand_prob = rand_prob,
    moderator = moderator,
    control = control,
    ml_method = "gam",
    cross_fit = FALSE,
    d_model_type = d_model,
    gam_control_spline_var = gam_control_spline_var
  )
  
  set.seed(1)
  
  fit_gam_cf <- emee_eif(
    dta = dta,
    id = id,
    outcome = outcome,
    treatment = treatment,
    rand_prob = rand_prob,
    moderator = moderator,
    control = control,
    ml_method = "gam",
    cross_fit = TRUE,
    d_model_type = d_model,
    gam_control_spline_var = gam_control_spline_var
  )
  
  set.seed(1)
  
  fit_rf_cf <- emee_eif(
    dta = dta,
    id = id,
    outcome = outcome,
    treatment = treatment,
    rand_prob = rand_prob,
    moderator = moderator,
    control = control,
    ml_method = "ranger",
    cross_fit = TRUE,
    d_model_type = d_model,
    gam_control_spline_var = gam_control_spline_var
  )
  
  set.seed(1)
  
  fit_stack_cf <- emee_eif(
    dta = dta,
    id = id,
    outcome = outcome,
    treatment = treatment,
    rand_prob = rand_prob,
    moderator = moderator,
    control = control,
    ml_method = "stack",
    cross_fit = TRUE,
    d_model_type = d_model,
    gam_control_spline_var = gam_control_spline_var
  )
  
  
  
  all_fits <- list(fit_emee,
                   fit_gam,
                   fit_gam_cf,
                   fit_rf_cf,
                   fit_stack_cf)
  names(all_fits) <- c("fit_emee",
                       "fit_gam",
                       "fit_gam_cf",
                       "fit_rf_cf",
                       "fit_stack_cf")
  
  dir.create("result", showWarnings = FALSE)
  filename <- paste0("count ", idesign, " - outcome = ", outcome, ", moderator = ", moderator)
  saveRDS(all_fits, file = paste0("result/", filename, ".RDS"))
}