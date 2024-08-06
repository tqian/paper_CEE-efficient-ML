eif_core <- function(
    dta,
    id,
    outcome,
    treatment,
    rand_prob,
    moderator, # vector of variables or NULL
    availability,
    mu_a0,
    mu_a1,
    d_vector = NULL, # 
    no_se = FALSE, # if TRUE, do not compute se (mostly for fitting beta_init)
    type = c("continuous", "binary")
    
) {
  sample_size <- length(unique(dta[[id]]))
  total_T <- nrow(dta) / sample_size
  
  Y <- dta[[outcome]]
  A <- dta[[treatment]]
  prob_A <- dta[[rand_prob]]
  S_mat <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator]))
  # S_mat is [S_{it}^T] stacked column by column, each column correspond to a (i,t) pair
  mu_hat_a0 <- dta[[mu_a0]]
  mu_hat_a1 <- dta[[mu_a1]]
  # Y_hat <- dta[[mu_a1]] * (1 - prob_A) + dta[[mu_a0]] * prob_A
  if (is.null(availability)) {
    avail <- rep(1, nrow(dta))
  } else {
    avail <- dta[[availability]]
  }
  
  if (is.null(d_vector)) {
    d_vector <- rep(1, total_T * sample_size)
  }
  
  Snames <- c("Intercept", moderator)
  dim_beta <- ncol(S_mat)
  
  # print(d_vector)
  
  # mu_hat_a0 <- ifelse(mu_hat_a0 ==0, 0.001, mu_hat_a0)
  # d_vector <- ifelse(d_vector >= 20, 20, d_vector)
  
  
  ee <- function(beta) {
    
    weight <- avail * (A - prob_A) / (prob_A * (1 - prob_A))
    if (type == "continuous") {
      residual <- Y - (A + prob_A - 1) * as.vector(S_mat %*% beta) - mu_hat_a1 * (1 - prob_A) - mu_hat_a0 * prob_A
    } else if (type == "binary") {
      residual <- exp(-1 * A * as.vector(S_mat %*% beta)) * Y - (1 - prob_A) * exp(- as.vector(S_mat %*% beta)) * mu_hat_a1 - 
        prob_A * mu_hat_a0
    }
    
    ef <- rep(NA, length(beta)) # value of estimating function
    for (ibeta in 1:dim_beta) {
      ef[ibeta] <- sum( d_vector * weight * residual * S_mat[, ibeta])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  solution <- tryCatch(
    {
      multiroot(ee, rep(0, dim_beta), useFortran = FALSE)
    },
    error = function(cond) {
      message(paste0("\nCatched error in multiroot with ee."))
      message(cond)
      return(list(root = rep(NaN, dim_beta), msg = cond,
                  f.root = rep(NaN, dim_beta)))
    })
  
  beta_hat <- solution$root
  
  if (no_se) {
    names(beta_hat) <- Snames
    return(list(beta_hat = beta_hat))
  } else {
    ## preparation for computing the variance estimator
    if (is.null(d_vector)) {
      d_vector <- rep(1, total_T * sample_size)
    }
    
    weight <- avail * (A - prob_A) / (prob_A * (1 - prob_A))
    
    if (type == "continuous") {
      residual <- Y - (A + prob_A - 1) * as.vector(S_mat %*% beta_hat) - mu_hat_a1 * (1 - prob_A) - mu_hat_a0 * prob_A
    } else if (type == "binary") {
      residual <- exp(-1 * A * as.vector(S_mat %*% beta_hat)) * Y - (1 - prob_A) * exp(- as.vector(S_mat %*% beta_hat)) * mu_hat_a1 - 
        prob_A * mu_hat_a0
    }
    
    
    
    
    total_person_decisionpoint <- total_T * sample_size
    r_term_collected <- residual
    D_term_collected <- matrix(NA, nrow = dim_beta, ncol = total_person_decisionpoint)
    partialr_partialbeta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = dim_beta)
    # print(d_vector)
    ## compute the meat term
    meat_term <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (i in 1:sample_size) {
      sum_dt_phit_for_this_i <- matrix(0, nrow = dim_beta, ncol = 1)
      for (t in 1:total_T) {
        it <- (i - 1) * total_T + t
        S_it <- matrix(S_mat[it, ], ncol = 1)
        dt_term <- d_vector[it]
        phit_term <- weight[it] * residual[it] * S_it
        if (type == "continuous") {
          partialr_partialbeta <- - (A[it] + prob_A[it] - 1) * S_it
        } else if (type == "binary") {
          partialr_partialbeta <- (- A[it] * exp(-1 * A[it] * as.vector(t(S_it) %*% beta_hat)) * Y[it] + 
            (1 - prob_A[it]) * exp(- as.vector(t(S_it) %*% beta_hat)) * mu_hat_a1[it]) * S_it
        }
        sum_dt_phit_for_this_i <- sum_dt_phit_for_this_i + dt_term * phit_term
        D_term_collected[, it] <- dt_term * weight[it] * S_it
        partialr_partialbeta_collected[it, ] <- partialr_partialbeta
      }
      meat_term <- meat_term + sum_dt_phit_for_this_i %*% t(sum_dt_phit_for_this_i)
    }
    meat_term <- meat_term / sample_size

    ## compute the bread term
    bread_term <- matrix(0, nrow = dim_beta, ncol = dim_beta)
    for (i in 1:sample_size) {
      for (t in 1:total_T) {
        it <- (i - 1) * total_T + t
        S_it <- matrix(S_mat[it, ], ncol = 1)
        dt_term <- d_vector[it]
        if (type == "continuous") {
          partial_phit_partial_beta_term <-
            - weight[it] * (A[it] + prob_A[it] - 1) * S_it %*% t(S_it)
        } else if (type == "binary") {
          partial_phit_partial_beta_term <-
            weight[it] * (-A[it] * exp(-1 * A[it] * as.vector(t(S_it) %*% beta_hat)) * Y[it] + 
            (1 - prob_A[it]) * exp(- as.vector(t(S_it) %*% beta_hat)) * mu_hat_a1[it]) * S_it %*% t(S_it)
        }
        # print(partial_phit_partial_beta_term)
        # print(mu_hat_a1[it])
        bread_term <- bread_term + dt_term * partial_phit_partial_beta_term
      }
    }
    bread_term <- bread_term / sample_size
    # print(bread_term)
    bread_term <- solve(bread_term, tol = 1e-50)
    
    beta_var <- bread_term %*% meat_term %*% t(bread_term) / sample_size
    beta_se <- sqrt(diag(beta_var))
    
    ## small sample correction (assume eta is fixed)
    
    person_first_index <- c(find_change_location(dta[, id]), total_person_decisionpoint + 1)
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
      D_term_i <- matrix(D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)], nrow = dim_beta)
      r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
      partialr_partialbeta_i <- partialr_partialbeta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
      H_ii <- partialr_partialbeta_i %*% bread_term %*% D_term_i / sample_size
      Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii, tol = 1e-50)
      
      Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- bread_term %*% Sigman_tilde %*% t(bread_term) / sample_size
    beta_se_adjusted <- sqrt(diag(varcov_adjusted))
    
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Snames
    
    ## calculate confidence interval
    
    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05/2, df = sample_size - dim_beta)
    conf_int_adjusted <- cbind(beta_hat - c * beta_se_adjusted,
                               beta_hat + c * beta_se_adjusted)
    colnames(conf_int) <- colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")
    
    
    return(list(beta_hat = beta_hat, 
                beta_se = beta_se,
                beta_se_adjusted = beta_se_adjusted,
                conf_int = conf_int,
                conf_int_adjusted = conf_int_adjusted))
  }
}

find_change_location <- function(v){
  n <- length(v)
  if (n <= 1) {
    stop("The vector need to have length > 1.")
  }
  return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}
