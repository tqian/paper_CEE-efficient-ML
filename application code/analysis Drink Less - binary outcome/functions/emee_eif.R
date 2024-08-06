# The one-step estimator
# (using estimating equation that does not require iterative estimation)

# Note 01.03:
# 1. Could refactor the code so that there is only on fit_d using the cf version and the
# non cf version could be implemented by setting dta_train and dta_holdout to be the same data set (whole data)
# 2. Could put the overlapping code in fit_d_part1 and fit_d_part2 in a function
# 3. could generalize the fit_d part to allow user input of residual and weight so that continuous and binary
# outcome could use the same function

# source("fit_d.R")
# source("eif.R")

stack.lm <- function(model_pred, # a data frame each column represents a prediction from weak learner
                     dta, 
                     outcome) {
  
  Y <- dta[[outcome]]
  X <- as.matrix(model_pred)
  nlearners <- ncol(model_pred)
  
  # glm and normalize the weights
  
  fit.glm <- glm(Y ~ X, family = binomial)
  weights <- fit.glm$coefficients
  weights[which(weights < 0 | is.na(weights))] <-  0
  weights <- as.vector(tail(weights, length(weights)-1))
  norm_weights <- weights / sum(weights)
  
  return(norm_weights)
  
}

get_predictions <- function(model_list, dta, col) {
  
  pred <- c()
  
  for (i in 1:length(model_list)) {
    regfit <- model_list[[i]]
    ml_method <- names(model_list)[i]
    colname <- paste0(col, ".", ml_method)
    
    if (ml_method == "ranger") {
      pred_fit <- list(as.vector(predict(regfit, data = dta)$predictions[, 2]))
      # pred_fit <- list(as.vector(predict(regfit, newdata = dta)))
    } else if (ml_method == "rf") {
      pred_fit <- list(as.vector(predict(regfit, newdata = dta, type = "prob")[, 2]))
    } else if (ml_method == "mean") {
      pred_fit <- list(as.vector(unlist(regfit)))
    } else {
      pred_fit <- list(as.vector(predict(regfit, newdata = dta, type = "response")))
    }
    
    names(pred_fit) <- colname
    pred <- c(pred, pred_fit)
  }
  
  return(data.frame(pred))
}

get_stack.model <- function(models, dta_a0, dta_a1, control, outcome, gam_control_spline_var) {
  
  stack.fit.a0 <- stack.fit.a1 <- c()
  
  for (i in 1:length(models)) {
    
    ml_method <- models[i]
    
    if (ml_method == "mean") {
      regfit_a0_mean <- list(mean = mean(as.numeric(dta_a0[[outcome]])))
      regfit_a1_mean <- list(mean = mean(as.numeric(dta_a1[[outcome]])))
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_mean)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_mean)
    } else if (ml_method == "gam") {
      
      if (length(gam_control_spline_var) == 0) {
        # gam_control_spline_var = NULL or gam_control_spline_var = c()
        gam_formula <- as.formula(paste0(outcome, " ~ ",
                                         paste0(control, collapse = " + ")))
      } else {
        gam_formula <- as.formula(
          paste0(outcome, " ~ ",
                 paste0(c(paste0("s(", gam_control_spline_var, ")"),
                          setdiff(control, gam_control_spline_var)),
                        collapse = " + ")
          )
        )
      }
      
      # gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
      regfit_a0 <- mgcv::gam(gam_formula, data = dta_a0, family = binomial)
      regfit_a1 <- mgcv::gam(gam_formula, data = dta_a1, family = binomial)
      # If we set binomial(link = "log"), it leads to error:
      #     Error in gam.fit3(x = X, y = y, sp = L %*% lsp1 + lsp0, Eb = Eb, UrS = UrS,  : 
      #     inner loop 2; can't correct step size
      regfit_a0_gam <- list(gam = regfit_a0)
      regfit_a1_gam <- list(gam = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_gam)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_gam)
    } else if (ml_method == "earth") {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- earth::earth(earth_formula, data = dta_a0, glm = list(family = binomial))
      regfit_a1 <- earth::earth(earth_formula, data = dta_a1, glm = list(family = binomial))
      regfit_a0_earth <- list(earth = regfit_a0)
      regfit_a1_earth <- list(earth = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_earth)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_earth)
    } else if (ml_method == "glm") {
      glm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- glm(glm_formula, data = dta_a0, family = binomial)
      regfit_a1 <- glm(glm_formula, data = dta_a1, family = binomial)
      regfit_a0_glm <- list(glm = regfit_a0)
      regfit_a1_glm <- list(glm = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_glm)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_glm)
    } else if (ml_method == "rf") {
      dta_a0[[outcome]] <- factor(dta_a0[[outcome]])
      dta_a1[[outcome]] <- factor(dta_a1[[outcome]])
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- randomForest::randomForest(rf_formula, data = dta_a0)
      regfit_a1 <- randomForest::randomForest(rf_formula, data = dta_a1)
      regfit_a0_rf <- list(rf = regfit_a0)
      regfit_a1_rf <- list(rf = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_rf)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_rf)
    } else if (ml_method == "ranger") {
      # Setting probability = TRUE to grow a probability forest as in Malley et al. (2012)
      dta_a0[[outcome]] <- factor(dta_a0[[outcome]])
      dta_a1[[outcome]] <- factor(dta_a1[[outcome]])
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- ranger::ranger(rf_formula, data = dta_a0, probability = TRUE)
      regfit_a1 <- ranger::ranger(rf_formula, data = dta_a1, probability = TRUE)
      
      # library(caret)
      # 
      # rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      # tune_grid <- expand.grid(mtry = seq(1, length(control), by = 1),
      #                          splitrule = "variance",
      #                          min.node.size = seq(1, 30, by = 5))
      # 
      # # A = 0
      # regfit_a0 <- train(rf_formula, method = "ranger", data = dta_a0,
      #                    tuneGrid = tune_grid,
      #                    trControl = trainControl(method = "cv", number = 5))
      # 
      # # A = 1
      # regfit_a1 <- train(rf_formula, method = "ranger", data = dta_a1,
      #                    tuneGrid = tune_grid,
      #                    trControl = trainControl(method = "cv", number = 5))
      
      regfit_a0_ranger <- list(ranger = regfit_a0)
      regfit_a1_ranger <- list(ranger = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_ranger)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_ranger)
    }
  }
  
  return(list(stack.fit.a0 = stack.fit.a0, stack.fit.a1 = stack.fit.a1))
  
}

fit_mu <- function(dta,
                   id,
                   outcome,
                   treatment,
                   rand_prob,
                   moderator, # vector of variables or NULL
                   control, # vector of variables
                   availability,
                   ml_method = c("gam", "earth", "earth_singlemodel",
                                 "earth_singlemodel2", "earth_singlemodel3",
                                 "lm",
                                 "rf", "ranger", "sl.smooth", "sl.all", 
                                 "stack"),
                   gam_control_spline_var
) {
  
  ml_method <- match.arg(ml_method)
  
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    "SL.ranger")
    # sl.library = c("SL.mean", "SL.gam", "SL.earth")
  } else if (ml_method == "stack") {
    stack.model <- c("gam", "earth", "mean")
  }
  
  dta_a0 <- dta[dta[[treatment]] == 0, ]
  dta_a1 <- dta[dta[[treatment]] == 1, ]
  
  if (ml_method == "gam") {
    
    if (length(gam_control_spline_var) == 0) {
      # gam_control_spline_var = NULL or gam_control_spline_var = c()
      gam_formula <- as.formula(paste0(outcome, " ~ ",
                                       paste0(control, collapse = " + ")))
    } else {
      gam_formula <- as.formula(
        paste0(outcome, " ~ ",
               paste0(c(paste0("s(", gam_control_spline_var, ")"),
                        setdiff(control, gam_control_spline_var)),
                      collapse = " + ")
        )
      )
    }
    
    # gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
    regfit_a0 <- mgcv::gam(gam_formula, data = dta_a0, family = binomial)
    regfit_a1 <- mgcv::gam(gam_formula, data = dta_a1, family = binomial)
    # If we set binomial(link = "log"), it leads to error:
    #     Error in gam.fit3(x = X, y = y, sp = L %*% lsp1 + lsp0, Eb = Eb, UrS = UrS,  : 
    #     inner loop 2; can't correct step size
  } else if (ml_method == "earth") {
    earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- earth::earth(earth_formula, data = dta_a0, glm = list(family = binomial))
    regfit_a1 <- earth::earth(earth_formula, data = dta_a1, glm = list(family = binomial))
  } else if (ml_method == "earth_singlemodel") {
    earth_formula <- as.formula(paste0(outcome, " ~ (", paste0(control, collapse = " + "), ") * ", treatment))
    regfit <- earth::earth(earth_formula, data = dta, glm = list(family = binomial))
  } else if (ml_method == "earth_singlemodel2") {
    earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ", treatment))
    regfit <- earth::earth(earth_formula, data = dta, glm = list(family = binomial))
  } else if (ml_method == "earth_singlemodel3") {
    if (is.null(moderator)) {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ", treatment))
    } else {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ",
                                         treatment, " * (", paste0(moderator, collapse = " + "), ")"))
    }
    regfit <- earth::earth(earth_formula, data = dta, glm = list(family = binomial))
  } else if (ml_method == "glm") {
    glm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- glm(glm_formula, data = dta_a0, family = binomial(link = "log"))
    regfit_a1 <- glm(glm_formula, data = dta_a1, family = binomial(link = "log"))
  } else if (ml_method == "rf") {
    dta_a0[[outcome]] <- factor(dta_a0[[outcome]])
    dta_a1[[outcome]] <- factor(dta_a1[[outcome]])
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- randomForest::randomForest(rf_formula, data = dta_a0)
    regfit_a1 <- randomForest::randomForest(rf_formula, data = dta_a1)
  } else if (ml_method == "ranger") {
    # Setting probability = TRUE to grow a probability forest as in Malley et al. (2012)
    dta_a0[[outcome]] <- factor(dta_a0[[outcome]])
    dta_a1[[outcome]] <- factor(dta_a1[[outcome]])
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- ranger::ranger(rf_formula, data = dta_a0, probability = TRUE)
    regfit_a1 <- ranger::ranger(rf_formula, data = dta_a1, probability = TRUE)
  } else if (ml_method %in% c("sl.smooth", "sl.all")) {
    regfit_a0 <- SuperLearner::SuperLearner(Y = dta_a0 %>% pull(!!outcome),
                                            X = dta_a0 %>% dplyr::select(!!control),
                                            family = binomial,
                                            verbose = FALSE,
                                            # id = dta_a0 %>% pull(!!id),
                                            SL.library = sl.library)
    regfit_a1 <- SuperLearner::SuperLearner(Y = dta_a1 %>% pull(!!outcome),
                                            X = dta_a1 %>% dplyr::select(!!control),
                                            family = binomial,
                                            verbose = FALSE,
                                            # id = dta_a1 %>% pull(!!id),
                                            SL.library = sl.library)
  } else if (ml_method == "stack") {
    
    # get weak learners
    fits <- get_stack.model(stack.model, 
                            dta_a0, 
                            dta_a1, 
                            control, 
                            outcome,
                            gam_control_spline_var)
    
    stack.fit.a0 <-fits$stack.fit.a0
    stack.fit.a1 <-fits$stack.fit.a1
    
    # get predictions from weak learners
    
    # A = 0
    pred_a0 <- get_predictions(stack.fit.a0, dta_a0, "eta_a0_hat")
    # print(pred_a0)
    
    # A = 1
    pred_a1 <- get_predictions(stack.fit.a1, dta_a1, "eta_a1_hat")
    
    # get weights
    weights_a0 <- stack.lm(pred_a0, dta_a0, outcome)
    # print(weights_a0)
    
    weights_a1 <- stack.lm(pred_a1, dta_a1, outcome)
    
    
    
  }
  
  # predicting eta_hat
  if (ml_method %in% c("sl.smooth", "sl.all")) {
    newdata_df <- dta %>% dplyr::select(!!control)
    dta <- dta %>%
      mutate(eta_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
             eta_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
  } else if (ml_method == "ranger") {
    dta <- dta %>%
      mutate(eta_hat_a0 = predict(regfit_a0, data = dta)$predictions[, 2],
             eta_hat_a1 = predict(regfit_a1, data = dta)$predictions[, 2])
  } else if (ml_method == "rf") {
    dta <- dta %>%
      mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta, type = "prob")[, 2],
             eta_hat_a1 = predict(regfit_a1, newdata = dta, type = "prob")[, 2])
  } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
    dta_seta0 <- dta_seta1 <- dta
    dta_seta0[[treatment]] <- 0
    dta_seta1[[treatment]] <- 1
    dta <- dta %>%
      mutate(eta_hat_a0 = predict(regfit, newdata = dta_seta0, type = "response"),
             eta_hat_a1 = predict(regfit, newdata = dta_seta1, type = "response"))
  } else if (ml_method == "stack"){
    
    # get predictions on data from each weak learner
    pred_a0 <- get_predictions(stack.fit.a0, dta, "eta_a0_hat")
    pred_a1 <- get_predictions(stack.fit.a1, dta, "eta_a1_hat")
    
    dta <- dta %>% 
      mutate(eta_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
             eta_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
    
  } else {
    dta <- dta %>%
      mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta, type = "response"),
             eta_hat_a1 = predict(regfit_a1, newdata = dta, type = "response"))
  }
  
  return(dta)
  
  
}

fit_mu_cf <- function(dta_train,
                      dta_holdout,
                      id,
                      outcome,
                      treatment,
                      rand_prob,
                      moderator, # vector of variables or NULL
                      control, # vector of variables
                      availability,
                      ml_method = c("gam", "earth", "earth_singlemodel",
                                    "earth_singlemodel2", "earth_singlemodel3",
                                    "lm",
                                    "rf", "ranger", "sl.smooth", "sl.all", 
                                    "stack"),
                      gam_control_spline_var
) {
  
  ml_method <- match.arg(ml_method)
  
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    "SL.ranger")
    # sl.library = c("SL.mean", "SL.gam", "SL.earth")
  } else if (ml_method == "stack") {
    stack.model <- c("gam", "earth", "mean")
  }
  
  # using cross-fitting as implemented in AIPW R package
  # (see their supplementary material for a illustrative diagram)
    
    dta_train_a0 <- dta_train[dta_train[[treatment]] == 0, ]
    dta_train_a1 <- dta_train[dta_train[[treatment]] == 1, ]
    
    if (ml_method == "gam") {
      
      if (length(gam_control_spline_var) == 0) {
        # gam_control_spline_var = NULL or gam_control_spline_var = c()
        gam_formula <- as.formula(paste0(outcome, " ~ ",
                                         paste0(control, collapse = " + ")))
      } else {
        gam_formula <- as.formula(
          paste0(outcome, " ~ ",
                 paste0(c(paste0("s(", gam_control_spline_var, ")"),
                          setdiff(control, gam_control_spline_var)),
                        collapse = " + ")
          )
        )
      }
      
      # gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
      regfit_a0 <- mgcv::gam(gam_formula, data = dta_train_a0, family = binomial)
      regfit_a1 <- mgcv::gam(gam_formula, data = dta_train_a1, family = binomial)
    } else if (ml_method == "earth") {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- earth::earth(earth_formula, data = dta_train_a0, glm = list(family = binomial))
      regfit_a1 <- earth::earth(earth_formula, data = dta_train_a1, glm = list(family = binomial))
    } else if (ml_method == "earth_singlemodel") {
      earth_formula <- as.formula(paste0(outcome, " ~ (", paste0(control, collapse = " + "), ") * ", treatment))
      regfit <- earth::earth(earth_formula, data = dta_train, glm = list(family = binomial))
    } else if (ml_method == "earth_singlemodel2") {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ", treatment))
      regfit <- earth::earth(earth_formula, data = dta_train, glm = list(family = binomial))
    } else if (ml_method == "earth_singlemodel3") {
      if (is.null(moderator)) {
        earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ", treatment))
      } else {
        earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + "), " + ",
                                           treatment, " * (", paste0(moderator, collapse = " + "), ")"))
      }
      regfit <- earth::earth(earth_formula, data = dta_train, glm = list(family = binomial))
    } else if (ml_method == "glm") {
      glm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- glm(glm_formula, data = dta_train_a0, family = binomial(link = "log"))
      regfit_a1 <- glm(glm_formula, data = dta_train_a1, family = binomial(link = "log"))
    } else if (ml_method == "rf") {
      dta_train_a0[[outcome]] <- factor(dta_train_a0[[outcome]])
      dta_train_a1[[outcome]] <- factor(dta_train_a1[[outcome]])
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- randomForest::randomForest(rf_formula, data = dta_train_a0)
      regfit_a1 <- randomForest::randomForest(rf_formula, data = dta_train_a1)
    } else if (ml_method == "ranger") {
      dta_train_a0[[outcome]] <- factor(dta_train_a0[[outcome]])
      dta_train_a1[[outcome]] <- factor(dta_train_a1[[outcome]])
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- ranger::ranger(rf_formula, data = dta_train_a0, probability = TRUE)
      regfit_a1 <- ranger::ranger(rf_formula, data = dta_train_a1, probability = TRUE)
    } else if (ml_method %in% c("sl.smooth", "sl.all")) {
      regfit_a0 <- SuperLearner::SuperLearner(Y = dta_train_a0 %>% pull(!!outcome),
                                              X = dta_train_a0 %>% dplyr::select(!!control),
                                              family = binomial,
                                              verbose = FALSE,
                                              id = dta_train_a0 %>% pull(!!id),
                                              SL.library = sl.library)
      regfit_a1 <- SuperLearner::SuperLearner(Y = dta_train_a1 %>% pull(!!outcome),
                                              X = dta_train_a1 %>% dplyr::select(!!control),
                                              family = binomial,
                                              verbose = FALSE,
                                              id = dta_train_a1 %>% pull(!!id),
                                              SL.library = sl.library)
    } else if (ml_method == "stack") {
      
      # get weak learners
      fits <- get_stack.model(stack.model, 
                              dta_train_a0, 
                              dta_train_a1, 
                              control, 
                              outcome,
                              gam_control_spline_var)
      
      stack.fit.a0 <-fits$stack.fit.a0
      stack.fit.a1 <-fits$stack.fit.a1
      
      # get predictions from weak learners
      
      # A = 0
      pred_a0 <- get_predictions(stack.fit.a0, dta_train_a0, "eta_a0_hat")
      # print(pred_a0)
      
      # A = 1
      pred_a1 <- get_predictions(stack.fit.a1, dta_train_a1, "eta_a1_hat")
      
      # get weights
      weights_a0 <- stack.lm(pred_a0, dta_train_a0, outcome)
      # print(weights_a0)
      
      weights_a1 <- stack.lm(pred_a1, dta_train_a1, outcome)
      
    }
    
    # predicting eta_hat for the holdout set
    if (ml_method %in% c("sl.smooth", "sl.all")) {
      newdata_df <- dta_holdout %>% dplyr::select(!!control)
      dta_holdout <- dta_holdout %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
               eta_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
    } else if (ml_method == "ranger") {
      dta_holdout <- dta_holdout %>%
        mutate(eta_hat_a0 = predict(regfit_a0, data = dta_holdout)$predictions[, 2],
               eta_hat_a1 = predict(regfit_a1, data = dta_holdout)$predictions[, 2])
    } else if (ml_method == "rf") {
      dta_holdout <- dta_holdout %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta_holdout, type = "prob")[, 2],
               eta_hat_a1 = predict(regfit_a1, newdata = dta_holdout, type = "prob")[, 2])
    } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
      dta_holdout_seta0 <- dta_holdout_seta1 <- dta_holdout
      dta_holdout_seta0[[treatment]] <- 0
      dta_holdout_seta1[[treatment]] <- 1
      dta_holdout <- dta_holdout %>%
        mutate(eta_hat_a0 = predict(regfit, newdata = dta_holdout_seta0, type = "response"),
               eta_hat_a1 = predict(regfit, newdata = dta_holdout_seta1, type = "response"))
    } else if (ml_method == "stack"){
      
      # get predictions on data from each weak learner
      pred_a0 <- get_predictions(stack.fit.a0, dta_holdout, "eta_a0_hat")
      pred_a1 <- get_predictions(stack.fit.a1, dta_holdout, "eta_a1_hat")
      
      dta_holdout <- dta_holdout %>% 
        mutate(eta_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
               eta_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
      
    } else {
      dta_holdout <- dta_holdout %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta_holdout, type = "response"),
               eta_hat_a1 = predict(regfit_a1, newdata = dta_holdout, type = "response"))
    }
    
    # predicting eta_hat for the train set
    # to use in training d term 
    if (ml_method %in% c("sl.smooth", "sl.all")) {
      newdata_df <- dta_train %>% dplyr::select(!!control)
      dta_train <- dta_train %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
               eta_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
    } else if (ml_method == "ranger") {
      dta_train <- dta_train %>%
        mutate(eta_hat_a0 = predict(regfit_a0, data = dta_train)$predictions[, 2],
               eta_hat_a1 = predict(regfit_a1, data = dta_train)$predictions[, 2])
    } else if (ml_method == "rf") {
      dta_train <- dta_train %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta_train, type = "prob")[, 2],
               eta_hat_a1 = predict(regfit_a1, newdata = dta_train, type = "prob")[, 2])
    } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
      dta_train_seta0 <- dta_train_seta1 <- dta_train
      dta_train_seta0[[treatment]] <- 0
      dta_train_seta1[[treatment]] <- 1
      dta_train <- dta_train %>%
        mutate(eta_hat_a0 = predict(regfit, newdata = dta_train_seta0, type = "response"),
               eta_hat_a1 = predict(regfit, newdata = dta_train_seta1, type = "response"))
    } else if (ml_method == "stack"){
      
      # get predictions on data from each weak learner
      pred_a0 <- get_predictions(stack.fit.a0, dta_train, "eta_a0_hat")
      pred_a1 <- get_predictions(stack.fit.a1, dta_train, "eta_a1_hat")
      
      dta_train <- dta_train %>% 
        mutate(eta_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
               eta_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
      
    } else {
      dta_train <- dta_train %>%
        mutate(eta_hat_a0 = predict(regfit_a0, newdata = dta_train, type = "response"),
               eta_hat_a1 = predict(regfit_a1, newdata = dta_train, type = "response"))
    }
    
    # dta_analysis <- rbind(dta_analysis, dta_holdout)
  
  return(list(dta_holdout = dta_holdout, dta_train = dta_train))
  
  
}

get_cf_fit <- function(id_folds, 
                       cf_fold,
                       id,
                       outcome,
                       treatment,
                       rand_prob,
                       moderator,
                       control,
                       availability,
                       ml_method,
                       d_model_type,
                       gam_control_spline_var) {
  
    for (k in 1:cf_fold) {

      # For each k, we want to get the following for dta_holdout:
      # 1. mu_k prediction on everyone in dta_holdout
      # 2. beta_init_k
      # 3. d_k prediction on everyone in dta_holdout
      #.   (this would depend on mu_k predictions and beta_init_k)

      ## Step 1: Sample split
      id_holdout <- na.omit(id_folds[, k])
      id_train <- na.omit(as.vector(id_folds[, -k]))

      dta_holdout <- dta[dta[[id]] %in% id_holdout, ]
      dta_train <- dta[dta[[id]] %in% id_train, ]

      ## Step 2: For each k, fit mu_k, beta_init_k, and d_k
      ##     Then obtain: prediction of mu_k on dta_holdout
      ##                  prediction on d_k on dta_holdout
      dta_train_and_holdout <- fit_mu_cf(
        dta_train,
        dta_holdout,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator,
        control,
        availability,
        ml_method,
        gam_control_spline_var
      )
      dta_train <- dta_train_and_holdout$dta_train
      dta_holdout <- dta_train_and_holdout$dta_holdout

      beta_init <- eif_core(
        dta_train,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator, # vector of variables or NULL
        availability,
        mu_a0 = "eta_hat_a0",
        mu_a1 = "eta_hat_a1",
        d_vector = NULL,
        no_se = TRUE,
        type = "binary"
      )$beta_hat

      d_array <- fit_d(
        dta_train,
        dta_holdout,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator,
        availability,
        mu_a0 = "eta_hat_a0",
        mu_a1 = "eta_hat_a1",
        beta = beta_init,
        d_model_type = d_model_type,
        outcome_type = "binary"
      )

      if (k == 1) {
        dta_collected <- dta_holdout
        d_array_collected <- d_array
        beta_init_collected <- matrix(beta_init, nrow = 1)
      } else {
        dta_collected <- rbind(dta_collected, dta_holdout)
        d_array_collected <- c(d_array_collected, d_array)
        beta_init_collected <- rbind(beta_init_collected, matrix(beta_init, nrow = 1))
      }
    }

    dta <- dta_collected
    d_array <- d_array_collected
    beta_init <- apply(beta_init_collected, 2, mean)
    
    fit <- eif_core(
      dta,
      id,
      outcome,
      treatment,
      rand_prob,
      moderator,
      availability,
      mu_a0 = "eta_hat_a0",
      mu_a1 = "eta_hat_a1",
      d_vector = d_array,
      no_se = FALSE,
      type = "binary"
    )
    
    return(fit)
    
}

emee_eif <- function(dta,
                     id,
                     outcome,
                     treatment,
                     rand_prob,
                     moderator, # vector of variables or NULL
                     control, # vector of variables
                     availability = NULL,
                     numerator_prob = 0.5,
                     ml_method = c("gam", "earth", "earth_singlemodel",
                                   "earth_singlemodel2", "earth_singlemodel3",
                                   "lm",
                                   "rf", "ranger", "sl.smooth", "sl.all",
                                   "stack"),
                     cross_fit = FALSE,
                     cf_fold = 10,
                     d_model_type = "empirical", # "empirical" - taking empirical average ignoring S or "linear-over-t" - pooled regression on S and t
                     return_d_array = FALSE, # if true, returns d_array for debugging purposes
                     user_specified_d = NULL, # specify a vector of length total_T if model_for_var = "user-specified-over-t"
                     return_beta_init = FALSE, # if true, returns beta_init
                     gam_control_spline_var = NULL
) {
    # gam: generalized additive model
    # glm: log-linear regression for binary outcome
    # rf: random forest (classification tree)
    # ranger: fast implementation of random forest (probability tree)
    # sl: super learner
    #   sl.smooth: sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
    #   sl.all:    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    #                             "SL.ranger", "SL.xgboost", "SL.nnet")
    
    ### Stage 1:
    ###     Fit mu_t, beta_init, d_t
    
    if (!cross_fit) {
      ### not cross-fitting
      
      ## Step 1: fit mu_t (this should be the same as the old code)
      #   Output: data set with columns (mu_1 and mu_0) attached
      dta <- fit_mu(dta,
                    id,
                    outcome,
                    treatment,
                    rand_prob,
                    moderator,
                    control,
                    availability,
                    ml_method,
                    gam_control_spline_var)
      
      ## Step 2: fit beta^init (by calling wcls_eif_core and ask not to calculate standard error)
      #   Output: beta_init
      beta_init <- eif_core(
        dta,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator, # vector of variables or NULL
        availability,
        mu_a0 = "eta_hat_a0",
        mu_a1 = "eta_hat_a1",
        d_vector = NULL,
        no_se = TRUE,
        type = "binary"
      )$beta_hat
      
      ## Step 3: fit d_t
      #   Step 3.1: fit the expected derivative
      #   Step 3.2: fit the expected variance
      d_array <- fit_d(
        dta,
        dta,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator,
        availability,
        mu_a0 = "eta_hat_a0",
        mu_a1 = "eta_hat_a1",
        beta = beta_init,
        d_model_type = d_model_type,
        outcome_type = "binary"
      )
      
      # print(d_array)
      # print(table(d_array))
      
    } else {
      ### cross-fitting
      
      # The following three lines of code (and the na.omit functions further down)
      # allows for situations where the number of individuals is not a multiple of cf_fold.
      
      # check whether there are both A=0 and A=1 in the cross fitted samples. if not, regenerate the samples.
      
      resplit <- TRUE
      
      while(resplit) {
      id_permuted <- sample(unique(dta[[id]]))
      id_permuted <- c(id_permuted, 
                       rep(NA, cf_fold - length(id_permuted) %% cf_fold))
      id_folds <- matrix(id_permuted, ncol = cf_fold, byrow = TRUE)
      
      out <- try(
        
        get_cf_fit(id_folds,
                   cf_fold,
                   id,
                   outcome,
                   treatment,
                   rand_prob,
                   moderator,
                   control,
                   availability,
                   ml_method,
                   d_model_type,
                   gam_control_spline_var) 

      )
      
      if (class(out) == "try-error") {
        resplit <- TRUE
      } else {
        resplit <- FALSE
        
        return(out)
      }
    }
}
      
      # count <- 0
      # for (i in 1:cf_fold) {
      #   dta_fold <- dta[dta[[id]] %in% na.omit(as.vector(id_folds[, -i])), ]
      #   
      #   dta_fold_a1 <- dta[dta[[treatment]] == 1, ]
      #   dta_fold_a0 <- dta[dta[[treatment]] == 0, ]
      #   
      #   if (nrow(dta_fold_a0) > 0 & nrow(dta_fold_a1) > 0) {
      #     count <- count + 1
      #   }
      #   
      # }
      
      # if (count == cf_fold) {
      # resplit <- FALSE
      # }
      
      
      
      # print(id_folds)
      
    #   for (k in 1:cf_fold) {
    #     
    #     # For each k, we want to get the following for dta_holdout:
    #     # 1. mu_k prediction on everyone in dta_holdout
    #     # 2. beta_init_k
    #     # 3. d_k prediction on everyone in dta_holdout
    #     #.   (this would depend on mu_k predictions and beta_init_k)
    #     
    #     ## Step 1: Sample split
    #     id_holdout <- na.omit(id_folds[, k])
    #     id_train <- na.omit(as.vector(id_folds[, -k]))
    #     
    #     dta_holdout <- dta[dta[[id]] %in% id_holdout, ]
    #     dta_train <- dta[dta[[id]] %in% id_train, ]
    #     
    #     ## Step 2: For each k, fit mu_k, beta_init_k, and d_k
    #     ##     Then obtain: prediction of mu_k on dta_holdout
    #     ##                  prediction on d_k on dta_holdout
    #     dta_train_and_holdout <- fit_mu_cf(
    #       dta_train,
    #       dta_holdout,
    #       id,
    #       outcome,
    #       treatment,
    #       rand_prob,
    #       moderator,
    #       control,
    #       availability,
    #       ml_method
    #     )
    #     dta_train <- dta_train_and_holdout$dta_train
    #     dta_holdout <- dta_train_and_holdout$dta_holdout
    #     
    #     beta_init <- eif_core(
    #       dta_train,
    #       id,
    #       outcome,
    #       treatment,
    #       rand_prob,
    #       moderator, # vector of variables or NULL
    #       availability,
    #       mu_a0 = "eta_hat_a0",
    #       mu_a1 = "eta_hat_a1",
    #       d_vector = NULL,
    #       no_se = TRUE,
    #       type = "binary"
    #     )$beta_hat
    #     
    #     d_array <- fit_d(
    #       dta_train,
    #       dta_holdout,
    #       id,
    #       outcome,
    #       treatment,
    #       rand_prob,
    #       moderator,
    #       availability,
    #       mu_a0 = "eta_hat_a0",
    #       mu_a1 = "eta_hat_a1",
    #       beta = beta_init,
    #       d_model_type = d_model_type,
    #       outcome_type = "binary"
    #     )
    #     
    #     if (k == 1) {
    #       dta_collected <- dta_holdout
    #       d_array_collected <- d_array
    #       beta_init_collected <- matrix(beta_init, nrow = 1)
    #     } else {
    #       dta_collected <- rbind(dta_collected, dta_holdout)
    #       d_array_collected <- c(d_array_collected, d_array)
    #       beta_init_collected <- rbind(beta_init_collected, matrix(beta_init, nrow = 1))
    #     }
    #   }
    #   
    #   dta <- dta_collected
    #   d_array <- d_array_collected
    #   beta_init <- apply(beta_init_collected, 2, mean)
    # }
    
    
    ### Stage 2: fit beta (by calling wcls_eif_core)
  
  # print(d_array)
    
    fit <- eif_core(
      dta,
      id,
      outcome,
      treatment,
      rand_prob,
      moderator,
      availability,
      mu_a0 = "eta_hat_a0",
      mu_a1 = "eta_hat_a1",
      d_vector = d_array,
      no_se = FALSE,
      type = "binary"
    )
    
    # beta_hat <- fit$beta_hat
    # beta_se <- fit$beta_se
    # Snames <- c("Intercept", moderator)
    # names(beta_hat) <- Snames
    # names(beta_se) <- Snames
    # names(beta_init) <- Snames
    # 
    # output <- list(beta_hat = beta_hat,
    #                beta_se = beta_se,
    #                beta_init = beta_init)
    
    if (return_beta_init) {
    output <- list(fit = fit, beta_init = beta_init)
    } else {
      output <- fit
    }

    if (return_d_array) {
      output <- c(output, list(d_array = d_array))
    }
    
    return(output)
}

emee_eif_oracle <- function(dta,
                       id,
                       outcome,
                       treatment,
                       rand_prob,
                       moderator, # vector of variables or NULL
                       control, # vector of variables
                       availability = NULL,
                       true_expect_y_a0,
                       true_expect_y_a1,
                       d_true
) {
  
  sample_size <- length(unique(dta[[id]]))
  total_T <- nrow(dta) / sample_size
  
  S_mat <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator]))
  dim_beta <- ncol(S_mat)
  
  # d_array_true <- array(NA, dim = c(nrow(dta), dim_beta, dim_beta))
  
  d_array_true <- matrix(NA, ncol = 1, nrow = nrow(dta))
  for (it in 1:nrow(dta)) {
    this_t <- (it - 1) %% total_T + 1
    # S_it <- matrix(S_mat[it, ], ncol = 1)
    if (length(d_true) == sample_size * total_T) {
      d_array_true[it, ] <- d_true[it]
    } else if (length(d_true) == total_T) {
      d_array_true[it, ] <- d_true[this_t]
    }
    
    # print(d_array_true)
    
    # The above line uses the same t-specific d_true value for all s
    # in E( phi_{it}^{otimes 2} | S_it = s ).
  }
    
    fit <- eif_core(dta = dta,
                    moderator = moderator,
                    treatment = treatment,
                    outcome = outcome,
                    mu_a0 = true_expect_y_a0,
                    mu_a1 = true_expect_y_a1,
                    rand_prob = rand_prob,
                    id = id,
                    availability = availability,
                    d_vector = d_array_true,
                    no_se = FALSE,
                    type = "binary")
    return(fit)
}

