# not centering at the residual
wcls_original <- function(dta,
                          id,
                          outcome,
                          treatment,
                          rand_prob,
                          moderator, # vector of variables or NULL
                          control, # vector of variables
                          availability = NULL,
                          numerator_prob = 0.5
) {
    
    if (is.character(numerator_prob)) {
        numerator_prob <- dta[[numerator_prob]]
    }
    
    outcome_var <- outcome
    control_var <- control
    trt_var <- treatment
    prob_A_var <- rand_prob
    id_var <- id
    prob_A_tilde <- numerator_prob
    
    prob_A <- dta[, prob_A_var]
    A <- dta[, trt_var]
    cA <- A - prob_A # centered A
    cA_tilde <- A - prob_A_tilde # centered A by \tilde{p}_t, used in WCLS
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator] ) )
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_var] ) )
    Y <- dta[, outcome_var]
    
    weight <- ifelse(A, prob_A_tilde / prob_A, 
                     (1 - prob_A_tilde) / (1 - prob_A))
    
    sample_size <- length(unique(dta[, id_var]))
    total_person_decisionpoint <- nrow(dta)
    
    Xnames <- c("Intercept", moderator)
    Znames <- c("Intercept", control_var)
    
    dim_alpha <- ncol(Zdm)
    dim_beta <- ncol(Xdm)
    
    p <- dim_beta
    q <- dim_alpha
    
    ee <- function(theta) {
        alpha <- as.matrix(theta[1:dim_alpha])
        beta <- as.matrix(theta[(dim_alpha + 1):(dim_alpha + dim_beta)])
        
        Zdm_alpha <- Zdm %*% alpha
        AXdm_beta <- A * (Xdm %*% beta)
        residual <- Y - AXdm_beta - Zdm_alpha
        weight <- ifelse(A, prob_A_tilde / prob_A, 
                         (1 - prob_A_tilde) / (1 - prob_A))
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:dim_alpha) {
            ef[i] <- sum( residual * weight * Zdm[, i])            
        }
        for (i in 1:dim_beta) {
            ef[i + dim_alpha] <- sum( residual * weight * cA_tilde * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    ##### Estimation #####
    
    solution <- tryCatch(
        {
            multiroot(ee, rep(0, dim_alpha + dim_beta), useFortran = FALSE)
        },
        error = function(cond) {
            message(paste0("\nCatched error in multiroot with ee."))
            message(cond)
            return(list(root = rep(NaN, dim_alpha + dim_beta), msg = cond,
                        f.root = rep(NaN, dim_alpha + dim_beta)))
        })
    alpha_hat <- solution$root[1:dim_alpha]
    beta_hat <- solution$root[(dim_alpha + 1):(dim_alpha + dim_beta)]
    
    ##### Asymptotic Variance #####
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        pre_multiplier <- weight[it]
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- 0
        partialD_partialtheta[1:q, (q+1):(q+p)] <- 0
        partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- 0
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - Zalpha - A[it] * Xbeta)
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- - c(Zdm[it, ], A[it] * Xdm[it, ])
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_var]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    ### 7. calculate confidence interval
    
    conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
    c <- qt(1 - 0.05/2, df = sample_size - p - q)
    conf_int_adjusted <- cbind(beta_hat - c * beta_se_adjusted,
                               beta_hat + c * beta_se_adjusted)
    colnames(conf_int) <- colnames(conf_int_adjusted) <- c("2.5 %", "97.5 %")
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                conf_int = conf_int, conf_int_adjusted = conf_int_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}


### code verifying that this is exactly the same as the MRTAnalysis package

if (0) {
    # first load a data set using simu_cont.R
    
    summary(MRTAnalysis::wcls(
        data = dta,
        id = "userid",
        outcome = "Y",
        treatment = "A",
        rand_prob = "prob_A",
        moderator_formula = ~ 1,
        control_formula = ~ S,
        availability = NULL,
        numerator_prob = 0.5,
        verbose = FALSE
    ))
    
    wcls_original(dta,
                  id = "userid",
                  outcome = "Y",
                  treatment = "A",
                  rand_prob = "prob_A",
                  moderator = NULL,
                  control = "S",
                  availability = NULL,
                  numerator_prob = 0.5
    )
}
