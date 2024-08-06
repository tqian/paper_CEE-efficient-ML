expit <- function(x){
  return(1/(1+exp(-x)))
}


const_rand_prob <- 0.5
const_S_halfrange <- 2
const_alpha0 <- 0.5
const_alpha3 <- -5
const_alpha1 <- 0.8
const_alpha2 <- 0.5

# const_alpha0 <- 2
# const_alpha1 <- 1
# const_alpha2 <- 0.5

const_beta0 <- 0.1
const_beta1 <- 0.2

# exp(-2.5) = 0.08
# exp(-1) = 0.37
# exp(-0.5) = 0.61
# exp(-0.2) = 0.82
# exp(-0.1) = 0.90


dgm_count_time_s_tvvar <- function(sample_size, total_T,
                              S_pattern = c("constant", "timevar"),
                              rand_prob_pattern = c("constant", "timevar"),
                              rand_prob_tuning_param = 1,
                              control_pattern = c("linear", "sine", "dbeta", "step"),
                              control_tuning_param = 0.5, # This should be between 0 and 1
                              ar1_tuning_param = 0.1 # how Y_{t+1} depends on Y_t; this should be between 0 and 0.1
) {
  rand_prob_pattern <- match.arg(rand_prob_pattern)
  control_pattern <- match.arg(control_pattern)
  
  df_names <- c("userid", "dp", "S",
                "A", "prob_A",
                "Y", "expect_Y", "expect_Y_A0", "expect_Y_A1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$dp <- rep(1:total_T, times = sample_size)
  
  A_lag1 <- rep(0, sample_size)
  prob_A_lag1 <- rep(0, sample_size)
  Y_lag1 <- rep(0, sample_size)
  
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
    
    # Everything covariate is transformed to within [0,1], so that
    # as long as control_tuning_param is well-controlled, we won't have 
    # success probability exceeding 1.
    
    # I will set the control part to be between exp(-2.5) and exp(-0.5)
    # exp(-2.5) = 0.08
    # exp(-0.5) = 0.61
    if (control_pattern == "linear") {
      S_transformed_within_unit_interval <- (S + const_S_halfrange + 1) / (2 * const_S_halfrange + 2)
      expect_Y_A0 <- 
        exp(
          const_alpha3 +
        const_alpha1 * dp +
        # const_alpha2 * S + 
          ar1_tuning_param * Y_lag1
      )
    } else if (control_pattern == "sine") {
      expect_Y_A0 <- 
        exp(
          const_alpha0 +
        control_pattern_tuning_param * sin(dp) +
        # control_pattern_tuning_param * sin(S) + 
          ar1_tuning_param * Y_lag1
      )
    } else if (control_pattern == "dbeta") {
      S_transformed_within_unit_interval <- (S + const_S_halfrange + 1) / (2 * const_S_halfrange + 2)
      expect_Y_A0 <- 
        exp(
          const_alpha0 +
        control_pattern_tuning_param * dbeta(dp / total_T, shape1 = 2, shape2 = 2) +
        # control_pattern_tuning_param * dbeta(S_transformed_within_unit_interval, shape1 = 2, shape2 = 2) +
        ar1_tuning_param * Y_lag1
      )
    } else if (control_pattern == "step") {
      expect_Y_A0 <- 
        exp(
          const_alpha0 +
        control_pattern_tuning_param * as.numeric(dp %% 2 == 0) +
        # control_pattern_tuning_param * as.numeric(round(S * 10) %% 2 == 0) +
        ar1_tuning_param * Y_lag1
      )
    }
    
    ## generate Y
    # expect_Y_A1 <- expect_Y_A0  + const_beta0 + const_beta1 * S
    expect_Y_A1 <- expect_Y_A0 * exp(const_beta0)
    expect_Y <- ifelse(A, expect_Y_A1, expect_Y_A0)
    Y <- rpois(sample_size, expect_Y)
    
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
    Y_lag1 <- Y
  }
  
  return(dta)
}