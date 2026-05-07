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
        ml_method = c("gam", "glm", "rf", "ranger", "sl.smooth", "sl.all", "sl.custom", "stack"),
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
      
      resplit <- TRUE
      
      while(resplit) {
        id_permuted <- sample(unique(dta[[id]]))
        id_permuted <- c(id_permuted, 
                         rep(NA, cf_fold - length(id_permuted) %% cf_fold))
        id_folds <- matrix(id_permuted, ncol = cf_fold, byrow = TRUE)
        
        out <- try(
          
          get_cf_fit(dta,
                     id_folds, 
                     cf_fold,
                     id,
                     outcome,
                     treatment,
                     rand_prob,
                     moderator,
                     control,
                     availability,
                     ml_method,
                     d_model_type) 
          
        )
        
        if (class(out) == "try-error") {
          resplit <- TRUE
        } else {
          resplit <- FALSE
          
          return(out)
        }
      }
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
    # 
    # if (return_d_array) {
    #     output <- c(output, list(d_array = d_array))
    # }
    
    return(fit)
}

get_cf_fit <- function(dta,
                       id_folds, 
                       cf_fold,
                       id,
                       outcome,
                       treatment,
                       rand_prob,
                       moderator,
                       control,
                       availability,
                       ml_method,
                       d_model_type) {
  
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
    
    # print(head(dta_train))
    # print(head(dta_holdout))
    
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
    
    d_array <- fit_d(
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
    mu_a0 = "mu_hat_a0",
    mu_a1 = "mu_hat_a1",
    d_vector = d_array,
    no_se = FALSE,
    type = "binary"
  )
  
  return(fit)
  
}

stack.lm <- function(model_pred, # a data frame each column represents a prediction from weak learner
                     dta, 
                     outcome) {
  
  Y <- dta[[outcome]]
  X <- as.matrix(model_pred)
  nlearners <- ncol(model_pred)
  
  # glm and normalize the weights
  
  fit.glm <- glm(Y ~ X, family = poisson)
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
      pred_fit <- list(as.vector(predict(regfit, data = dta, type = "response")$predictions))
      
    } else {
      pred_fit <- list(as.vector(predict(regfit, newdata = dta, type = "response")))
      
    }
    names(pred_fit) <- colname
    pred <- c(pred, pred_fit)
  }
  
  return(data.frame(pred))
}

get_stack.model <- function(models, dta_a0, dta_a1, control, outcome) {
  
  stack.fit.a0 <- stack.fit.a1 <- c()
  
  for (i in 1:length(models)) {
    
    ml_method <- models[i]
    
    if (ml_method == "gam") {
      gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
      regfit_a0 <- gam(gam_formula, data = dta_a0, family = poisson)
      regfit_a1 <- gam(gam_formula, data = dta_a1, family = poisson)
      regfit_a0_gam <- list(gam = regfit_a0)
      regfit_a1_gam <- list(gam = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_gam)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_gam)
    } else if (ml_method == "glm") {
      lm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- glm(lm_formula, data = dta_a0, family = poisson)
      regfit_a1 <- glm(lm_formula, data = dta_a1, family = poisson)
      regfit_a0_glm <- list(glm = regfit_a0)
      regfit_a1_glm <- list(glm = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_glm)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_glm)
    } else if (ml_method == "earth") {
      earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- earth::earth(earth_formula, data = dta_a0, glm = list(family = poisson))
      regfit_a1 <- earth::earth(earth_formula, data = dta_a1, glm = list(family = poisson))
      regfit_a0_earth <- list(earth = regfit_a0)
      regfit_a1_earth <- list(earth = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_earth)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_earth)
    } else if (ml_method == "rf") {
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- randomForest(rf_formula, data = dta_a0)
      regfit_a1 <- randomForest(rf_formula, data = dta_a1)
      regfit_a0_rf <- list(rf = regfit_a0)
      regfit_a1_rf <- list(rf = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_rf)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_rf)
    } else if (ml_method == "ranger") {
      rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
      regfit_a0 <- ranger(rf_formula, data = dta_a0)
      regfit_a1 <- ranger(rf_formula, data = dta_a1)
      regfit_a0_ranger <- list(ranger = regfit_a0)
      regfit_a1_ranger <- list(ranger = regfit_a1)
      stack.fit.a0 <- c(stack.fit.a0, regfit_a0_ranger)
      stack.fit.a1 <- c(stack.fit.a1, regfit_a1_ranger)
    }
  }
  
  return(list(stack.fit.a0 = stack.fit.a0, stack.fit.a1 = stack.fit.a1))
  
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
    ml_method = c("gam", "glm", "rf", "ranger", "sl.smooth", "sl.all", "sl.custom", "stack")
) {
  
  # gam: generalized additive model
  # lm: linear regression
  # ranger: fast implementation of random forest
  # sl: super learner
  
  ml_method <- match.arg(ml_method)
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
                   "SL.ranger")
  } else if (ml_method == "sl.custom") {
    
    SL.glm.poisson <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
    {
      if (is.matrix(X)) {
        X = as.data.frame(X)
      }
      fit.glm <- glm(Y ~ ., data = X, family = poisson, weights = obsWeights, 
                     model = model)
      if (is.matrix(newX)) {
        newX = as.data.frame(newX)
      }
      pred <- predict(fit.glm, newdata = newX, type = "response")
      fit <- list(object = fit.glm)
      class(fit) <- "SL.glm"
      out <- list(pred = pred, fit = fit)
      return(out)
    }
    
    SL.gam.poisson <- function (Y, X, newX, family, obsWeights, cts.num = 4, 
                                ...) 
    {
      if ("mgcv" %in% loadedNamespaces()) 
        warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
      cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
      if (sum(!cts.x) > 0) {
        gam.model <- as.formula(paste0("Y~", paste0(paste0("s(", 
                                                        colnames(X[, cts.x, drop = FALSE]), 
                                                        ")"), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop = FALSE]), collapse = "+")))
      }
      else {
        gam.model <- as.formula(paste0("Y~", paste0(paste0("s(", 
                                                        colnames(X[, cts.x, drop = FALSE]), 
                                                        ")"), collapse = "+")))
      }
      if (sum(!cts.x) == length(cts.x)) {
        gam.model <- as.formula(paste0("Y~", paste0(colnames(X), 
                                                  collapse = "+")))
      }
      
      fit.gam <- gam::gam(gam.model, data = X, family = poisson, 
                          control = gam::gam.control(maxit = 50, bf.maxit = 50), 
                          weights = obsWeights)
      if (packageVersion("gam") >= "1.15") {
        pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
      }
      else {
        stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
      }

      fit <- list(object = fit.gam)
      out <- list(pred = pred, fit = fit)
      class(out$fit) <- c("SL.gam")
      return(out)
    }
    
    sl.library = c("SL.mean", "SL.glm.poisson", "SL.gam.poisson", "SL.earth",
                   "SL.ranger")
  } else if (ml_method == "stack") {
    stack.model <- c("gam", "earth", "glm", "ranger")
  }
  
  dta_a0 <- dta[dta[[treatment]] == 0, ]
  dta_a1 <- dta[dta[[treatment]] == 1, ]
  
  if (ml_method == "gam") {
    gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
    regfit_a0 <- gam(gam_formula, data = dta_a0, family = poisson)
    regfit_a1 <- gam(gam_formula, data = dta_a1, family = poisson)
  } else if (ml_method == "glm") {
    lm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- glm(lm_formula, data = dta_a0, family = poisson)
    regfit_a1 <- glm(lm_formula, data = dta_a1, family = poisson)
  } else if (ml_method == "earth") {
    earth_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- earth::earth(earth_formula, data = dta_a0, glm = list(family = poisson))
    regfit_a1 <- earth::earth(earth_formula, data = dta_a1, glm = list(family = poisson))
  } else if (ml_method == "rf") {
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- randomForest(rf_formula, data = dta_a0)
    regfit_a1 <- randomForest(rf_formula, data = dta_a1)
  } else if (ml_method == "ranger") {
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- ranger(rf_formula, data = dta_a0)
    regfit_a1 <- ranger(rf_formula, data = dta_a1)
  } else if (ml_method %in% c("sl.smooth", "sl.all", "sl.custom")) {
    
    regfit_a0 <- SuperLearner(Y = dta_a0 %>% pull(!!outcome),
                              X = dta_a0 %>% dplyr::select(!!control),
                              family = gaussian(),
                              verbose = FALSE,
                              SL.library = sl.library)
    regfit_a1 <- SuperLearner(Y = dta_a1 %>% pull(!!outcome),
                              X = dta_a1 %>% dplyr::select(!!control),
                              family = gaussian(),
                              verbose = FALSE,
                              SL.library = sl.library)
  } else if (ml_method == "stack") {
    
    # get weak learners
    fits <- get_stack.model(stack.model, 
                            dta_a0, 
                            dta_a1, 
                            control, 
                            outcome)
    
    stack.fit.a0 <-fits$stack.fit.a0
    stack.fit.a1 <-fits$stack.fit.a1
    
    # get predictions from weak learners
    
    # A = 0
    pred_a0 <- get_predictions(stack.fit.a0, dta_a0, "mu_a0_hat")
    # print(pred_a0)
    
    # A = 1
    pred_a1 <- get_predictions(stack.fit.a1, dta_a1, "mu_a1_hat")
    
    # get weights
    weights_a0 <- stack.lm(pred_a0, dta_a0, outcome)
    # print(weights_a0)
    
    weights_a1 <- stack.lm(pred_a1, dta_a1, outcome)
    
    
    
  }
  
  # predicting mu_hat
  if (ml_method %in% c("sl.smooth", "sl.all", "sl.custom")) {
    newdata_df <- dta %>% dplyr::select(!!control)
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, type = "response")$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, type = "response")$pred)
  } else if (ml_method == "ranger") {
    dta <- dta %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta, type = "response")$predictions,
             mu_hat_a1 = predict(regfit_a1, data = dta, type = "response")$predictions)
  } else if (ml_method == "stack"){
    
    # get predictions on data from each weak learner
    pred_a0 <- get_predictions(stack.fit.a0, dta, "mu_a0_hat")
    pred_a1 <- get_predictions(stack.fit.a1, dta, "mu_a1_hat")
    
    dta <- dta %>% 
      mutate(mu_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
             mu_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
    
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
    ml_method = c("gam", "glm", "rf", "ranger", "sl.smooth", "sl.all", "sl.custom", "stack")
) {
  
  # gam: generalized additive model
  # lm: linear regression
  # ranger: fast implementation of random forest
  # sl: super learner
  
  ml_method <- match.arg(ml_method)
  if (ml_method == "sl.smooth") {
    sl.library <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  } else if (ml_method == "sl.all") {
    sl.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth",
                   "SL.ranger")
  } else if (ml_method == "sl.custom") {
    
    SL.glm.poisson <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
    {
      if (is.matrix(X)) {
        X = as.data.frame(X)
      }
      fit.glm <- glm(Y ~ ., data = X, family = poisson, weights = obsWeights, 
                     model = model)
      if (is.matrix(newX)) {
        newX = as.data.frame(newX)
      }
      pred <- predict(fit.glm, newdata = newX, type = "response")
      fit <- list(object = fit.glm)
      class(fit) <- "SL.glm"
      out <- list(pred = pred, fit = fit)
      return(out)
    }
    
    SL.glm.poisson <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) 
    {
      if (is.matrix(X)) {
        X = as.data.frame(X)
      }
      fit.glm <- glm(Y ~ ., data = X, family = poisson, weights = obsWeights, 
                     model = model)
      if (is.matrix(newX)) {
        newX = as.data.frame(newX)
      }
      pred <- predict(fit.glm, newdata = newX, type = "response")
      fit <- list(object = fit.glm)
      class(fit) <- "SL.glm"
      out <- list(pred = pred, fit = fit)
      return(out)
    }
    
    SL.gam.poisson <- function (Y, X, newX, family, obsWeights, cts.num = 4, 
                                ...) 
    {
      if ("mgcv" %in% loadedNamespaces()) 
        warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
      cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
      if (sum(!cts.x) > 0) {
        gam.model <- as.formula(paste0("Y~", paste0(paste0("s(", 
                                                           colnames(X[, cts.x, drop = FALSE]), 
                                                           ")"), collapse = "+"), "+", paste(colnames(X[, !cts.x, drop = FALSE]), collapse = "+")))
      }
      else {
        gam.model <- as.formula(paste0("Y~", paste0(paste0("s(", 
                                                           colnames(X[, cts.x, drop = FALSE]), 
                                                           ")"), collapse = "+")))
      }
      if (sum(!cts.x) == length(cts.x)) {
        gam.model <- as.formula(paste0("Y~", paste0(colnames(X), 
                                                    collapse = "+")))
      }
      
      fit.gam <- gam::gam(gam.model, data = X, family = poisson, 
                          control = gam::gam.control(maxit = 50, bf.maxit = 50), 
                          weights = obsWeights)
      if (packageVersion("gam") >= "1.15") {
        pred <- gam::predict.Gam(fit.gam, newdata = newX, type = "response")
      }
      else {
        stop("This SL.gam wrapper requires gam version >= 1.15, please update the gam package with 'update.packages('gam')'")
      }
      
      fit <- list(object = fit.gam)
      out <- list(pred = pred, fit = fit)
      class(out$fit) <- c("SL.gam")
      return(out)
    }
    
    sl.library = c("SL.mean", "SL.glm.poisson", "SL.gam.poisson", "SL.earth",
                   "SL.ranger")
  } else if (ml_method == "stack") {
    stack.model <- c("gam", "earth", "glm", "ranger")
  }
  
  ## Fit mu_k using B_k^c (dta_train)
  dta_train_a0 <- dta_train[dta_train[[treatment]] == 0, ]
  dta_train_a1 <- dta_train[dta_train[[treatment]] == 1, ]
  
  if (ml_method == "gam") {
    gam_formula <- as.formula(paste0(outcome, " ~ ", paste0("s(", control, ")", collapse = " + ")))
    regfit_a0 <- gam(gam_formula, data = dta_train_a0, family = poisson)
    regfit_a1 <- gam(gam_formula, data = dta_train_a1, family = poisson)
  } else if (ml_method == "glm") {
    lm_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- glm(lm_formula, data = dta_train_a0, family = poisson)
    regfit_a1 <- glm(lm_formula, data = dta_train_a1, family = poisson)
  } else if (ml_method == "rf") {
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- randomForest(rf_formula, data = dta_train_a0)
    regfit_a1 <- randomForest(rf_formula, data = dta_train_a1)
  } else if (ml_method == "ranger") {
    rf_formula <- as.formula(paste0(outcome, " ~ ", paste0(control, collapse = " + ")))
    regfit_a0 <- ranger(rf_formula, data = dta_train_a0)
    regfit_a1 <- ranger(rf_formula, data = dta_train_a1)
  } else if (ml_method %in% c("sl.smooth", "sl.all", "sl.custom")) {
    regfit_a0 <- SuperLearner(Y = dta_train_a0 %>% pull(!!outcome),
                              X = dta_train_a0 %>% dplyr::select(!!control),
                              family = gaussian(),
                              verbose = FALSE,
                              SL.library = sl.library)
    regfit_a1 <- SuperLearner(Y = dta_train_a1 %>% pull(!!outcome),
                              X = dta_train_a1 %>% dplyr::select(!!control),
                              family = gaussian(),
                              verbose = FALSE,
                              SL.library = sl.library)
    
  } else if (ml_method == "stack") {
    
    # get weak learners
    fits <- get_stack.model(stack.model, 
                            dta_train_a0, 
                            dta_train_a1, 
                            control, 
                            outcome)
    
    stack.fit.a0 <-fits$stack.fit.a0
    stack.fit.a1 <-fits$stack.fit.a1
    
    # get predictions from weak learners
    
    # A = 0
    pred_a0 <- get_predictions(stack.fit.a0, dta_train_a0, "mu_a0_hat")
    # print(pred_a0)
    
    # A = 1
    pred_a1 <- get_predictions(stack.fit.a1, dta_train_a1, "mu_a1_hat")
    
    # get weights
    weights_a0 <- stack.lm(pred_a0, dta_train_a0, outcome)
    # print(weights_a0)
    
    weights_a1 <- stack.lm(pred_a1, dta_train_a1, outcome)
    
  }
  
  ## Predict mu_hat on B_k (dta_holdout)
  ## This will be used in the final estimating equation
  if (ml_method %in% c("sl.smooth", "sl.all", "sl.custom")) {
    newdata_df <- dta_holdout %>% dplyr::select(!!control)
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, onlySL = T)$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, onlySL = T)$pred)
    
  } else if (ml_method == "ranger") {
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta_holdout, type = "response")$predictions,
             mu_hat_a1 = predict(regfit_a1, data = dta_holdout, type = "response")$predictions)
  } else if (ml_method == "stack"){
    
    # get predictions on data from each weak learner
    pred_a0 <- get_predictions(stack.fit.a0, dta_holdout, "mu_a0_hat")
    pred_a1 <- get_predictions(stack.fit.a1, dta_holdout, "mu_a1_hat")
    
    dta_holdout <- dta_holdout %>% 
      mutate(mu_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
             mu_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
    
  } else {
    dta_holdout <- dta_holdout %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = dta_holdout, type = "response"),
             mu_hat_a1 = predict(regfit_a1, newdata = dta_holdout, type = "response"))
  }
  
  ## Predict mu_hat on B_k^c (dta_train)
  ## This will be used in computing d_array later on
  if (ml_method %in% c("sl.smooth", "sl.all", "sl.custom")) {
    newdata_df <- dta_train %>% dplyr::select(!!control)
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, newdata = newdata_df, onlySL = T)$pred,
             mu_hat_a1 = predict(regfit_a1, newdata = newdata_df, onlySL = T)$pred)

  } else if (ml_method == "ranger") {
    dta_train <- dta_train %>%
      mutate(mu_hat_a0 = predict(regfit_a0, data = dta_train, type = "response")$predictions,
             mu_hat_a1 = predict(regfit_a1, data = dta_train, type = "response")$predictions)
  } else if (ml_method == "stack"){
    
    # get predictions on data from each weak learner
    pred_a0 <- get_predictions(stack.fit.a0, dta_train, "mu_a0_hat")
    pred_a1 <- get_predictions(stack.fit.a1, dta_train, "mu_a1_hat")
    
    dta_train <- dta_train %>% 
      mutate(mu_hat_a0 = as.matrix(pred_a0) %*% weights_a0,
             mu_hat_a1 = as.matrix(pred_a1) %*% weights_a1)
    
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