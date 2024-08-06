library(mvtnorm)

expit <- function(x){
  return(1/(1+exp(-x)))
}

exponential_corr <- function(base, total_T) {
    vcov <- matrix(NA, nrow = total_T, ncol = total_T)
    for (i in 1:nrow(vcov)) {
        for (j in 1:nrow(vcov)) {
            vcov[i, j] <- base^(abs(j - i) / 2)
        }
    }
    return(vcov)
}

const_beta0 <- 0.5
const_beta1 <- 0.2

const_rand_prob <- 0.5
const_S_halfrange <- 2
const_alpha0 <- 1
const_alpha1 <- 1
const_alpha2 <- 1


dgm_cont_time_s_tvvar <- function(sample_size,
                            total_T,
                            S_pattern = c("constant", "timevar"),
                            rand_prob_pattern = c("constant", "timevar"),
                            rand_prob_tuning_param = 1,
                            control_pattern = c("linear", "sine", "dbeta", "step"),
                            control_pattern_tuning_param = 1,
                            error_var_pattern = c("ind", "timevar-linear"), 
                            error_var_pattern_tuning_param = 1, # controls how non-constant the variance is over time
                            error_var_scale_tuning_param = 1, # variance at time 1
                            error_corr_tuning_param = 0 # base of the exponentially-decaying correlation, set to 0 for independent error terms
) {
    rand_prob_pattern <- match.arg(rand_prob_pattern)
    control_pattern <- match.arg(control_pattern)
    
    df_names <- c("userid", "dp", "S",
                  "A", "prob_A",
                  "Y", "expect_Y", "expect_Y_A0", "expect_Y_A1", "eps")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$dp <- rep(1:total_T, times = sample_size)
    
    error_corr_mat <- exponential_corr(error_corr_tuning_param, total_T)
    
    if (error_var_pattern == "ind") {
        error_sd_vec <- rep(error_var_scale_tuning_param, total_T)
    } else if (error_var_pattern == "timevar-linear") {
        error_sd_vec <- seq(from = error_var_scale_tuning_param, 
                             by = error_var_pattern_tuning_param, length = total_T)
    }
    error_sd_vec <- sqrt(error_sd_vec)
    
    error_sd_diag <- diag(c(error_sd_vec))
    error_vcov_mat <- error_sd_diag %*% error_corr_mat %*% error_sd_diag
    
    dta$eps <- as.vector(t(rmvnorm(n = sample_size,
                                   mean = rep(0, nrow(error_vcov_mat)),
                                   sigma = error_vcov_mat)))
    # TQ checked on 2023.11.20: the above code does the correct vectorization
    # (concatenate within each individual first, then across individuals)
    
    A_lag1 <- rep(0, sample_size)
    prob_A_lag1 <- rep(0, sample_size)
    
    for (t in 1:total_T) {
        # row index for the rows corresponding to dp t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        
        dp <- dta$dp[row_index]
        
        ## generate S
        if (S_pattern == "constant") {
            S <- runif(sample_size, min = - const_S_halfrange, max = const_S_halfrange)
        } else if (S_pattern == "timevar") {
            S <- runif(sample_size, min = - const_S_halfrange, max = const_S_halfrange) + A_lag1 - prob_A_lag1
        }
        
        ## generate A
        if (rand_prob_pattern == "constant") {
            prob_A <- rep(const_rand_prob, sample_size)
        } else if (rand_prob_pattern == "timevar") {
            S_transformed <- S / (2 * const_S_halfrange + 2)
            dp_transformed <- (dp - total_T/2) / total_T
            prob_A <- expit(rand_prob_tuning_param * S_transformed + 
                                rand_prob_tuning_param * dp_transformed)
            prob_A <- pmax(pmin(prob_A, 0.9), 0.1)
        }
        A <- rbinom(sample_size, 1, prob = prob_A)
        
        ## generate control part
        if (control_pattern == "linear") {
            expect_Y_A0 <- 
                const_alpha0 + 
                const_alpha1 * dp + 
                const_alpha2 * S
        } else if (control_pattern == "sine") {
            expect_Y_A0 <- 
                const_alpha0 + 
                control_pattern_tuning_param * sin(dp) + 
                control_pattern_tuning_param * sin(S)
        } else if (control_pattern == "dbeta") {
            S_transformed_within_unit_interval <- (S + const_S_halfrange + 1) / (2 * const_S_halfrange + 2)
            expect_Y_A0 <- 
                const_alpha0 + 
                control_pattern_tuning_param * dbeta(dp / total_T, shape1 = 2, shape2 = 2) + 
                control_pattern_tuning_param * dbeta(S_transformed_within_unit_interval, shape1 = 2, shape2 = 2)
        } else if (control_pattern == "step") {
            expect_Y_A0 <- 
                const_alpha0 + 
                control_pattern_tuning_param * as.numeric(dp %% 2 == 0) +
                control_pattern_tuning_param * as.numeric(round(S * 10) %% 2 == 0)
        }
        
        ## generate Y
        expect_Y_A1 <- expect_Y_A0 + const_beta0 + const_beta1 * S
        expect_Y <- ifelse(A, expect_Y_A1, expect_Y_A0)
        
        Y <- expect_Y + dta$eps[row_index]
        
        ## Put variables into dta
        dta$S[row_index] <- S
        dta$A[row_index] <- A
        dta$prob_A[row_index] <- prob_A
        dta$Y[row_index] <- Y
        dta$expect_Y[row_index] <- expect_Y
        dta$expect_Y_A0[row_index] <- expect_Y_A0
        dta$expect_Y_A1[row_index] <- expect_Y_A1
        
        # Gather lagged variables to be used
        A_lag1 <- A
        prob_A_lag1 <- prob_A
    }
    
    return(dta)
}