# functions for fitting d part

train_d_part2 <- function(dta_train,
                                   id,
                                   outcome,
                                   treatment,
                                   rand_prob,
                                   moderator,
                                   availability,
                                   mu_a0,
                                   mu_a1,
                                   beta, # usually this is beta_init
                                   d_model_type, # "empirical" or "linear-over-t"
                         outcome_type # "continuous" or "binary"
) {
  sample_size <- length(unique(dta_train[[id]]))
  total_T <- nrow(dta_train) / sample_size
  
  # mu_a0 = "expect_Y_A0"
  # mu_a1 = "expect_Y_A1"
  # beta = 0.5
  
  Y <- dta_train[[outcome]]
  A <- dta_train[[treatment]]
  prob_A <- dta_train[[rand_prob]]
  S_mat <- as.matrix(cbind( rep(1, nrow(dta_train)), dta_train[, moderator]))
  # S_mat is [S_{it}^T] stacked column by column, each column correspond to a (i,t) pair
  mu_hat_a0 <- dta_train[[mu_a0]]
  mu_hat_a1 <- dta_train[[mu_a1]]
  if (is.null(availability)) {
    avail <- rep(1, nrow(dta_train))
  } else {
    avail <- dta_train[[availability]]
  }
  
  dim_beta <- ncol(S_mat)
  
  if (outcome_type == "binary") {
    ### 1. Compute the value of phi_{it}^{otimes 2}
    
    residual <- exp(-1 * A * as.vector(S_mat %*% beta)) * Y - (1 - prob_A) * exp(- as.vector(S_mat %*% beta)) * mu_hat_a1 - 
      prob_A * mu_hat_a0
    weight <- avail * (A - prob_A) / (prob_A * (1 - prob_A))
    
    phi_otimes2 <- residual^2 * weight^2
  } else if (outcome_type == "continuous") {
    ### 1. Compute the value of phi_{it}^{otimes 2}
    
    Y_hat <- dta_train[[mu_a1]] * (1 - prob_A) + dta_train[[mu_a0]] * prob_A
    
    residual <- Y - (A + prob_A - 1) * as.vector(S_mat %*% beta) - Y_hat
    weight <- avail * (A - prob_A) / (prob_A * (1 - prob_A))
    
    phi_otimes2 <- residual^2 * weight^2
  }

  ### 2. Estimate E[phi_{it}^{otimes 2} | S_it]
  if (d_model_type == "linear-over-t") {
    predictor <- cbind(dta_train[, moderator], rep(1:total_T, times = sample_size), phi_otimes2)
    pred_df <- as.data.frame(predictor)
    colnames(pred_df) <- c(moderator, "t", "phi_otimes2")
    fit_formula <- as.formula(paste0("phi_otimes2~", paste0(moderator, collapse = "+"), "+t"))
    fit <- lm(fit_formula, data = pred_df)
    
    return(fit)
  } else if (d_model_type == "empirical") {
    ## t-specific estimate
    t_specific_estimate <- matrix(NA, nrow = total_T, ncol = 1)
    
    for (t in 1:total_T) {
      it_idx <-  seq(from = t, by = total_T, length = sample_size)
      t_specific_estimate[t, ] <- mean(phi_otimes2[it_idx])
    }
    
    return(t_specific_estimate)
  } else if (d_model_type == "earth") {
    predictor <- cbind(dta_train[, moderator], rep(1:total_T, times = sample_size), phi_otimes2)
    pred_df <- as.data.frame(predictor)
    colnames(pred_df) <- c(moderator, "t", "phi_otimes2")
    # if (is.null(moderator) | sum(dta_train[, moderator] != rep(1:total_T, times = sample_size)) == 0) {
    #   fit_formula <- as.formula(paste0("phi_otimes2~", "t"))
    # } else {
      fit_formula <- as.formula(paste0("phi_otimes2~", paste0(moderator, collapse = "+"), "+t"))
    # }
    fit <- earth::earth(fit_formula, data = pred_df, degree = 10)
    # print(summary(fit))
    return(fit)
  } else if (d_model_type == "gam") {
    predictor <- cbind(dta_train[, moderator], rep(1:total_T, times = sample_size), phi_otimes2)
    pred_df <- as.data.frame(predictor)
    colnames(pred_df) <- c(moderator, "t", "phi_otimes2")
    if (is.null(moderator)) {
      fit_formula <- as.formula(paste0("phi_otimes2~", "+s(t)"))
    } else {
      fit_formula <- as.formula(paste0("phi_otimes2~", paste0("s(", moderator, ")", collapse = "+"), "+s(t)"))
    }
    fit <- mgcv::gam(fit_formula, data = pred_df)
    return(fit)
  } else if (d_model_type == "rf") {
    predictor <- cbind(dta_train[, moderator], rep(1:total_T, times = sample_size), phi_otimes2)
    pred_df <- as.data.frame(predictor)
    colnames(pred_df) <- c(moderator, "t", "phi_otimes2")
    fit_formula <- as.formula(paste0("phi_otimes2~", paste0(moderator, collapse = "+"), "+t"))
    fit <- randomForest::randomForest(fit_formula, data = pred_df)
    return(fit)
  }
  
}

predict.dpart1cts <- function(object, newdata) {
  data_rows <- nrow(newdata)
  return(rep(-1, data_rows))
}

fit_d_part_all <- function(dta_holdout,
                        id,
                        outcome,
                        treatment,
                        rand_prob,
                        moderator,
                        availability,
                        mu_a0,
                        mu_a1,
                        beta, # usually this is beta_init
                        fit_train, # the fitted model from dta_train
                        part, # 1 or 2 corresponding to the first and second term in d array -- combining overlapping code for debugging purpose
                        d_model_type # "empirical" or "linear-over-t"
) {
  
  sample_size_holdout <- length(unique(dta_holdout[[id]]))
  total_T <- nrow(dta_holdout) / sample_size_holdout
  S_mat_holdout <- as.matrix(cbind( rep(1, nrow(dta_holdout)), dta_holdout[, moderator]))
  
  if (d_model_type %in% c("linear-over-t", "earth", "gam", "rf")) {
    # predict on holdout set for linear-over-t model
    predictor_holdout <- cbind(dta_holdout[, moderator], rep(1:total_T, times = sample_size_holdout))
    predictor_holdout <- as.data.frame(predictor_holdout)
    colnames(predictor_holdout) <- c(moderator, "t")
    pred_values <- predict(fit_train, newdata = predictor_holdout)
    
    # pred_values_holdout <- matrix(NA, nrow = nrow(dta_holdout), ncol = 1)
    # 
    # for (t in 1:total_T) {
    #   it_idx <-  seq(from = t, by = total_T, length = sample_size_holdout)
    #   pred_values_holdout[it_idx, ] <- mean(pred_values[it_idx])
    # }
    
    pred_values_holdout <- pred_values
    
  
  } else if (d_model_type == "empirical") {
    # predict on holdout set for taking empirical average
    if (length(fit_train) == 1) {
      pred_values_holdout <- rep(fit_train, times = nrow(dta_holdout))
    } else if (length(fit_train) == total_T) {
      pred_values_holdout <- matrix(NA, nrow = nrow(dta_holdout), ncol = 1)
      for (t in 1:total_T) {
        it_idx <- seq(from = t, by = total_T, length = sample_size_holdout)
        pred_values_holdout[it_idx, ] <- fit_train[t]
      }
    } else if (length(fit_train) == total_T * sample_size_holdout) {
      pred_values_holdout <- fit_train
    }
    
  } else {
    stop("Please input a valid d model.")
  }
  
  if (part == 1) {
    d_part <- as.vector(pred_values_holdout)
  } else if (part == 2) {
  d_part <- 1 / as.vector(pred_values_holdout)
  } else {
    stop("Please enter either 1 or 2 to indicate first or second term of d array accordingly.")
  }
  
  return(d_part)
}

train_d_part1 <- function(
    dta_train, 
    id,
    outcome,
    treatment,
    rand_prob,
    moderator,
    availability,
    mu_a0,
    mu_a1,
    beta,
    outcome_type,
    d_model_type
) {
  if (outcome_type == "continuous") {
    if (d_model_type == "empirical") {
      return(-1)
    } else if (d_model_type %in% c("linear-over-t", "earth", "gam", "rf")) {
      fit <- function() {
        return(-1)
      }
      class(fit) <- "dpart1cts"
      return(fit)
    }
    # return(-1)
  } else if (outcome_type == "binary") {
    sample_size <- length(unique(dta_train[[id]]))
    total_T <- nrow(dta_train) / sample_size
    
    # mu_a0 = "expect_Y_A0"
    # mu_a1 = "expect_Y_A1"
    # beta = 0.5
    
    Y <- dta_train[[outcome]]
    A <- dta_train[[treatment]]
    prob_A <- dta_train[[rand_prob]]
    S_mat <- as.matrix(cbind( rep(1, nrow(dta_train)), dta_train[, moderator]))
    # S_mat is [S_{it}^T] stacked column by column, each column correspond to a (i,t) pair
    mu_hat_a0 <- dta_train[[mu_a0]]
    mu_hat_a1 <- dta_train[[mu_a1]]
    if (is.null(availability)) {
      avail <- rep(1, nrow(dta_train))
    } else {
      avail <- dta_train[[availability]]
    }
    
    dim_beta <- ncol(S_mat)
    
    ### 1. Compute the value of phi_{it}^{otimes 2}
    
    residual <- -A * exp(-1 * A * as.vector(S_mat %*% beta)) * Y + (1 - prob_A) * exp(- as.vector(S_mat %*% beta)) * mu_hat_a1
    weight <- avail * (A - prob_A) / (prob_A * (1 - prob_A))
    
    phi <- residual * weight
    
    ### 2. Estimate E[phi_{it}^{otimes 2} | S_it]
    if (d_model_type == "linear-over-t") {
      predictor <- cbind(S_mat, rep(1:total_T, times = sample_size), phi)
      pred_df <- as.data.frame(predictor)
      colnames(pred_df) <- c(moderator, "t", "phi")
      fit_formula <- as.formula(paste0("phi~", paste0(moderator, collapse = "+"), "+t"))
      fit <- lm(fit_formula, data = pred_df)
      # fit <- lm(phi ~ predictor)
      return(fit)
    } else if (d_model_type == "empirical") {
      return(mean(phi))
    }
    
    # return(mean(phi))
  } else {
    stop("Please input a valid outcome type.")
  }
}

fit_d <- function(
    dta_train, # train set
    dta_holdout, # holdout set, for non cross fit, dta_train = dta_holdout = dta
    id,
    outcome,
    treatment,
    rand_prob,
    moderator,
    availability,
    mu_a0,
    mu_a1,
    beta, # usually this is beta_init
    d_model_type, # either "empirical" or "linear-over-t"
    outcome_type # either "continuous" or "binary"
) {

    ### train d_array part 1
  
  fit_part1 <- train_d_part1(
    dta_train, 
    id,
    outcome,
    treatment,
    rand_prob,
    moderator,
    availability,
    mu_a0,
    mu_a1,
    beta,
    outcome_type,
    d_model_type
  )
  
  ## predict d_array part 1 on holdout set
    
    d_part1 <- fit_d_part_all(
      dta_holdout,
      id,
      outcome,
      treatment,
      rand_prob,
      moderator,
      availability,
      mu_a0,
      mu_a1,
      beta,
      fit_part1,
      part = 1,
      d_model_type
    )
    
  ## train d_array part 2
    
    fit_part2 <- train_d_part2(
      dta_train, 
      id,
      outcome,
      treatment,
      rand_prob,
      moderator,
      availability,
      mu_a0,
      mu_a1,
      beta,
      d_model_type,
      outcome_type
    )
    
    ## predict d_array part 2 on holdout set
    
    d_part2 <- fit_d_part_all(
      dta_holdout,
      id,
      outcome,
      treatment,
      rand_prob,
      moderator,
      availability,
      mu_a0,
      mu_a1,
      beta,
      fit_part2,
      part = 2,
      d_model_type
    )
    
    d_array <- d_part1 * d_part2
  
  return(d_array)
}
