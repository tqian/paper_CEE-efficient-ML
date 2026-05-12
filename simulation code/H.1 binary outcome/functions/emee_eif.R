# Note 2023.11.19: 
# 1. currently does not incorporate nontrivial availability
#    (e.g., in fit_mu, did not subset the data set to remove unavailable decision points)
# 2. currently does not allow specification of "no intercept" in either moderator or control
#    (e.g., in fit_mu, did not allow adding "-1" in the formula)
#    (e.g., in wcls_eif_core, S_mat is constructed by always adding the intercept)
# 3. currently each individual must have equal number of decision points
#    (e.g., in fit_d_part2, total_T is computed by nrow(dta) / sample_size)
# 4. currently the fitting of d is super simplistic:
#.   for the derivative, uses a single number (avg_weight)
#.   for the expected variance, ignores S and computes a empirical average for each t separately
# 5. currently requires dta to be sorted by id then decision point

# We could simplify the d term to - 1 / E(res_t^2 w_t^2 | S_t)

# source("fit_d.R")
# source("eif.R")

# The EIF-based estimator
emee_eif <- function(
        dta,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator, # vector of variables or NULL
        control, # vector of variables
        availability = NULL,
        ml_method = c("gam", "earth", "earth_singlemodel",
                      "earth_singlemodel2", "earth_singlemodel3",
                      "lm",
                      "rf", "ranger", "sl.smooth", "sl.all"),
        cross_fit = FALSE,
        cf_fold = 10,
        d_model_type = "empirical", # "empirical" - taking empirical average ignoring S or "linear-over-t" - pooled regression on S and t
        return_d_array = FALSE, # if true, returns d_array for debugging purposes
        user_specified_d = NULL # specify a vector of length total_T if model_for_var = "user-specified-over-t"
) {
    
    # model_for_var <- match.arg(model_for_var)
    ml_method <- match.arg(ml_method)
    
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
                      ml_method)
        
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
            mu_a0 = "mu_hat_a0",
            mu_a1 = "mu_hat_a1",
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
            mu_a0 = "mu_hat_a0",
            mu_a1 = "mu_hat_a1",
            beta = beta_init,
            d_model_type = d_model_type,
            outcome_type = "binary"
        )
        
        # print(table(d_array))
        
    } else {
      ### cross-fitting
      
      # The following three lines of code (and the na.omit functions further down)
      # allows for situations where the number of individuals is not a multiple of cf_fold.
      id_permuted <- sample(unique(dta[[id]]))
      id_permuted <- c(id_permuted, 
                       rep(NA, cf_fold - length(id_permuted) %% cf_fold))
      id_folds <- matrix(id_permuted, ncol = cf_fold, byrow = TRUE)
      
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
          ml_method
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
          mu_a0 = "mu_hat_a0",
          mu_a1 = "mu_hat_a1",
          d_vector = NULL,
          no_se = TRUE,
          type = "binary"
        )$beta_hat
        
        d_array_k <- fit_d(
          dta_train,
          dta_holdout,
          id,
          outcome,
          treatment,
          rand_prob,
          moderator,
          availability,
          mu_a0 = "mu_hat_a0",
          mu_a1 = "mu_hat_a1",
          beta = beta_init,
          d_model_type = d_model_type,
          outcome_type = "binary"
        )
        
        
        if (k == 1) {
          dta_collected <- dta_holdout
          d_array_collected <- d_array_k
          beta_init_collected <- matrix(beta_init, nrow = 1)
        } else {
          dta_collected <- rbind(dta_collected, dta_holdout)
          d_array_collected <- c(d_array_collected, d_array_k)
          beta_init_collected <- rbind(beta_init_collected, matrix(beta_init, nrow = 1))
        }
      }
      
      dta <- dta_collected
      d_array <- d_array_collected
      beta_init <- apply(beta_init_collected, 2, mean)
    }
    
    # print(d_array)

    
    ### Stage 2: fit beta (by calling wcls_eif_core)
    
    fit <- eif_core(
        dta,
        id,
        outcome,
        treatment,
        rand_prob,
        moderator,
        availability,
        mu_a0 = "mu_hat_a0",
        mu_a1 = "mu_hat_a1",
        d_vector = d_array,
        no_se = FALSE,
        type = "binary"
    )
    
    return(fit)
}


fit_mu <- function(
    dta,
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
                  "rf", "ranger", "sl.smooth", "sl.all")
) {
  
  # gam: generalized additive model
  # lm: linear regression
  # ranger: fast implementation of random forest
  # sl: super learner
  
  ml_method <- match.arg(ml_method)
  
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    # sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    # "SL.ranger", "SL.xgboost", "SL.nnet")
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
                   "SL.ranger")
  }
  
  dta_a0 <- dta[dta[[treatment]] == 0, ]
  dta_a1 <- dta[dta[[treatment]] == 1, ]
  
  if (ml_method == "gam") {
    
    gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
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
  }
  
  # predicting mu_hat
  if (ml_method %in% c("sl.smooth", "sl.all")) {
    newdata_df <- dta %>% dplyr::select(!!control)
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
  } else if (ml_method == "ranger") {
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta)$predictions[, 2],
             mu_hat_a1 = predict(regfit_a1, data = dta)$predictions[, 2])
  } else if (ml_method == "rf") {
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta, type = "prob")[, 2],
             mu_hat_a1 = predict(regfit_a1, newdata = dta, type = "prob")[, 2])
  } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
    dta_seta0 <- dta_seta1 <- dta
    dta_seta0[[treatment]] <- 0
    dta_seta1[[treatment]] <- 1
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit, newdata = dta_seta0, type = "response"),
             mu_hat_a1 = predict(regfit, newdata = dta_seta1, type = "response"))
  } else {
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta, type = "response"),
             mu_hat_a1 = predict(regfit_a1, newdata = dta, type = "response"))
  }
  
  return(dta)
}


fit_mu_cf <- function(
    dta_train, # B_k^c
    dta_holdout, # B_k
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
                  "rf", "ranger", "sl.smooth", "sl.all")
) {
  
  # gam: generalized additive model
  # lm: linear regression
  # ranger: fast implementation of random forest
  # sl: super learner
  
  ml_method <- match.arg(ml_method)
  
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    # sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
    # "SL.ranger", "SL.xgboost", "SL.nnet")
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
                   "SL.ranger")
  }
  
  # using cross-fitting as implemented in AIPW R package
  # (see their supplementary material for a illustrative diagram)
  
  dta_train_a0 <- dta_train[dta_train[[treatment]] == 0, ]
  dta_train_a1 <- dta_train[dta_train[[treatment]] == 1, ]
  
  if (ml_method == "gam") {
    
    gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
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
  }
  
  # predicting mu_hat for the holdout set
  if (ml_method %in% c("sl.smooth", "sl.all")) {
    newdata_df <- dta_holdout %>% dplyr::select(!!control)
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
  } else if (ml_method == "ranger") {
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta_holdout)$predictions[, 2],
             mu_hat_a1 = predict(regfit_a1, data = dta_holdout)$predictions[, 2])
  } else if (ml_method == "rf") {
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta_holdout, type = "prob")[, 2],
             mu_hat_a1 = predict(regfit_a1, newdata = dta_holdout, type = "prob")[, 2])
  } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
    dta_holdout_seta0 <- dta_holdout_seta1 <- dta_holdout
    dta_holdout_seta0[[treatment]] <- 0
    dta_holdout_seta1[[treatment]] <- 1
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit, newdata = dta_holdout_seta0, type = "response"),
             mu_hat_a1 = predict(regfit, newdata = dta_holdout_seta1, type = "response"))
  } else {
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta_holdout, type = "response"),
             mu_hat_a1 = predict(regfit_a1, newdata = dta_holdout, type = "response"))
  }
  
  # predicting mu_hat for the train set
  # to use in training d term 
  if (ml_method %in% c("sl.smooth", "sl.all")) {
    newdata_df <- dta_train %>% dplyr::select(!!control)
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
  } else if (ml_method == "ranger") {
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta_train)$predictions[, 2],
             mu_hat_a1 = predict(regfit_a1, data = dta_train)$predictions[, 2])
  } else if (ml_method == "rf") {
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta_train, type = "prob")[, 2],
             mu_hat_a1 = predict(regfit_a1, newdata = dta_train, type = "prob")[, 2])
  } else if (ml_method %in% c("earth_singlemodel", "earth_singlemodel2", "earth_singlemodel3")) {
    dta_train_seta0 <- dta_train_seta1 <- dta_train
    dta_train_seta0[[treatment]] <- 0
    dta_train_seta1[[treatment]] <- 1
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit, newdata = dta_train_seta0, type = "response"),
             mu_hat_a1 = predict(regfit, newdata = dta_train_seta1, type = "response"))
  } else {
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta_train, type = "response"),
             mu_hat_a1 = predict(regfit_a1, newdata = dta_train, type = "response"))
  }
  
  return(list(dta_train = dta_train,
              dta_holdout = dta_holdout))
}

emee_eif_oracle <- function(dta,
                            id,
                            outcome,
                            treatment,
                            rand_prob,
                            moderator, # vector of variables or NULL
                            control, # vector of variables
                            availability = NULL,
                            mu_a0_true,
                            mu_a1_true,
                            d_true # a vector of length total T (only handles the most simple case of S = \emptyset)
) {
  sample_size <- length(unique(dta[[id]]))
  total_T <- nrow(dta) / sample_size
  
  S_mat <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator]))
  dim_beta <- ncol(S_mat)
  
  d_array_true <- matrix(NA, ncol = 1, nrow = nrow(dta))
  for (it in 1:nrow(dta)) {
    this_t <- (it - 1) %% total_T + 1
    # S_it <- matrix(S_mat[it, ], ncol = 1)
    if (length(d_true) == sample_size * total_T) {
      d_array_true[it, ] <- d_true[it]
    } else if (length(d_true) == total_T) {
      d_array_true[it, ] <- d_true[this_t]
    } else if (length(d_true) == 1) {
      d_array_true[it, ] <- d_true
    }
  }
  
  fit <- eif_core(
    dta,
    id,
    outcome,
    treatment,
    rand_prob,
    moderator,
    availability,
    mu_a0 = mu_a0_true,
    mu_a1 = mu_a1_true,
    d_vector = d_array_true,
    no_se = FALSE,
    type = "binary"
  )
  
  return(fit)
}

# get_cf_fit <- function(dta,
#                        id_folds, 
#                        cf_fold,
#                        id,
#                        outcome,
#                        treatment,
#                        rand_prob,
#                        moderator,
#                        control,
#                        availability,
#                        ml_method,
#                        d_model_type) {
#   
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
#     # print(head(dta_train))
#     # print(head(dta_holdout))
#     
#     beta_init <- eif_core(
#       dta_train,
#       id,
#       outcome,
#       treatment,
#       rand_prob,
#       moderator, # vector of variables or NULL
#       availability,
#       mu_a0 = "mu_hat_a0",
#       mu_a1 = "mu_hat_a1",
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
#       mu_a0 = "mu_hat_a0",
#       mu_a1 = "mu_hat_a1",
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
#   
#   fit <- eif_core(
#     dta,
#     id,
#     outcome,
#     treatment,
#     rand_prob,
#     moderator,
#     availability,
#     mu_a0 = "mu_hat_a0",
#     mu_a1 = "mu_hat_a1",
#     d_vector = d_array,
#     no_se = FALSE,
#     type = "binary"
#   )
#   
#   return(fit)
#   
# }