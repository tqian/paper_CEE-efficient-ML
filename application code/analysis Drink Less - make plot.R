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
mypalette <- brewer.pal(n = 7, "Set1")

design <- expand.grid(moderator = c("NULL", "decision_index"),
                           outcome = c("seconds_8to9", "count_8to9", "binary_8to9"),
                           stringsAsFactors = FALSE)

set.seed(1)

design$idesign <- ifelse(design$moderator == "NULL", 1, 3)
design$outcome_type <- case_when(design$outcome == "seconds_8to9" ~ "cont",
                                 design$outcome == "count_8to9" ~ "count",
                                 design$outcome == "binary_8to9" ~ "bin")

# plot_titles <- c("Outcome: continuous \n Moderator: None",
#                  "Outcome: continuous \n Moderator: decision index")

plot_est_collected <- c()
plot_re_collected <- c()
plot_collected <- c()

for (i in 1:nrow(design)) {
  idesign <- design$idesign[i]
  outcome <- design$outcome[i]
  moderator <- design$moderator[i]
  if (moderator == "NULL") {
    moderator <- NULL
  }
  outcome_type <- design$outcome_type[i]
  filename <- paste0(outcome_type, " ", idesign, " - outcome = ", outcome, ", moderator = ", moderator)
  foldername <- paste0("analysis Drink Less - ", case_when(outcome_type == "cont" ~ "continuous", 
                          outcome_type == "count" ~ "count", 
                          outcome_type == "bin" ~ "binary"), " outcome")
  
  all_fits <- readRDS(paste0(foldername, "result/", filename, ".RDS"))
  outcome_title <- case_when(outcome_type == "cont" ~ "continuous",
                             outcome_type == "bin" ~ "binary",
                             outcome_type == "count" ~ "count")
  moderator_title <- ifelse(is.null(moderator), "None", "decision index")
  plot_title <- paste0("Outcome: ", outcome_title, "\n Moderator: ", moderator_title)
  
  # WCLS, SL.CF, GAM.CF, GAM, RF.CF
  # if (outcome_type == "count") {
  #   all_fits <- all_fits[c(1, 5, 2:4)]
  #   print(all_fits)
  # }
  # all_fits <- all_fits[c(1, 2, 3, 5, 7, 8, 9)] # only keep the non-cf ones because we haven't implemented the cf for the new eif
  
  
  est <- t(sapply(all_fits, function(x) x$beta_hat))
  se <- t(sapply(all_fits, function(x) x$beta_se_adjusted))
  
  if (is.null(moderator)) {
    
    re <- t((t(se[2:length(se)]) / se[1])^(-2))
    est1 <- as.numeric(est)
    se1 <- as.numeric(se)
    lci1 <- est1 - 1.96 * se1
    rci1 <- est1 + 1.96 * se1
    re1 <- as.numeric(re)
    result <- tibble(method = as.factor(names(all_fits)),
                     est = est1, se = se1, lci = lci1, rci = rci1,
                     re = c(1, re1))
    
    relabel_function <- function(x) {
      y <- case_when(
        x == "fit_wcls" ~ "WCLS",
        x == "fit_emee" ~ "EMEE",
        x == "fit_sl_cf" ~ "SL.CF",
        x == "fit_sl" ~ "SL",
        x == "fit_gam_cf" ~ "GAM.CF",
        x == "fit_gam" ~ "GAM",
        x == "fit_rf_cf" ~ "RF.CF",
        x == "fit_rf" ~ "RF",
        x == "fit_stack_cf" ~ "STACK.CF",
        x == "fit_stack" ~ "STACK",
        TRUE ~ "unknown"
      )
      return(y)
    }
    
    result <- result %>%
      mutate(method = fct_relabel(method, relabel_function)) %>%
      mutate(method = fct_relevel(method,
                                  ifelse(outcome_type == "cont", "WCLS", "EMEE"),
                                  "SL.CF", "SL",
                                  "STACK.CF", "STACK",
                                  "GAM.CF", "GAM",
                                  "RF.CF", "RF"))
    
    print(result)
    print(filename)
    
    ylab_label_p1 <- TeX(r'($\hat{\beta}_0$ ($\pm$1.96 SE))')
    
    ylab_label_p2 <- TeX(r'(RE($\hat{\beta}_0$))')
    
    p1 <- result %>%
      ggplot(aes(x = method, y = est, color = method)) + 
      geom_errorbar(aes(ymin = lci, ymax = rci), width = .2, size = 1.2) +
      geom_point(size = 4) +
      theme_bw() +
      scale_color_manual(guide = "none", values = mypalette) +
      xlab(NULL) +
      ylab(ylab_label_p1) +
      theme(axis.text.x = element_text(size = 11, color = "black", face = "bold"),
            axis.text.y = element_text(size = 11),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12, hjust = 0.5))
    p1 <- p1 + ggtitle(plot_title)
    
    p2 <- result %>%
      # filter(!(method %in% c("SL", "RF"))) %>% 
      ggplot(aes(x = method, y = re, group = 1)) +
      geom_point(aes(color = method), size = 4) +
      # geom_line() +
      geom_hline(yintercept = 1, linetype = 2) +
      coord_cartesian(ylim = c(0.8, 2.3)) +
      #scale_y_continuous(n.breaks = 3) +
      theme_bw() +
      scale_color_manual(guide = "none", values = mypalette) +
      xlab(NULL) +
      ylab(ylab_label_p2) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 11))
    
    # blank_plot <- ggplot() + 
    #   theme(panel.background = element_rect(fill = 'white', color = 'white'))
    
    plot_est_collected <- c(plot_est_collected, list(p1))
    plot_re_collected <- c(plot_re_collected, list(p2))
    # plot_est_collected <- c(plot_est_collected, list(blank_plot))
    # plot_re_collected <- c(plot_re_collected, list(blank_plot))
    
  } else {
    re <- t((t(se[2:nrow(se), ]) / se[1, ])^(-2))
    
    for (k in 1:2) {
      est1 <- est[, k]
      se1 <- se[, k]
      lci1 <- est1 - 1.96 * se1
      rci1 <- est1 + 1.96 * se1
      re1 <- re[, k]
      result <- tibble(method = as.factor(names(all_fits)),
                       est = est1, se = se1, lci = lci1, rci = rci1,
                       re = c(1, re1))
      
      relabel_function <- function(x) {
        y <- case_when(
          x == "fit_wcls" ~ "WCLS",
          x == "fit_emee" ~ "EMEE",
          x == "fit_sl_cf" ~ "SL.CF",
          x == "fit_sl" ~ "SL",
          x == "fit_gam_cf" ~ "GAM.CF",
          x == "fit_gam" ~ "GAM",
          x == "fit_rf_cf" ~ "RF.CF",
          x == "fit_rf" ~ "RF",
          x == "fit_stack_cf" ~ "STACK.CF",
          x == "fit_stack" ~ "STACK",
          TRUE ~ "unknown"
        )
        return(y)
      }
      print(result)
      print(filename)
      result <- result %>% 
        mutate(method = fct_relabel(method, relabel_function)) %>% 
        mutate(method = fct_relevel(method,
                                    ifelse(outcome_type == "cont", "WCLS", "EMEE"),
                                    "SL.CF", "SL",
                                    "STACK.CF", "STACK",
                                    "GAM.CF", "GAM",
                                    "RF.CF", "RF"))
      
      ylab_label_p1 <- ifelse(k == 1,
                              TeX(r'($\hat{\beta}_0$ ($\pm$1.96 SE))'),
                              TeX(r'($\hat{\beta}_1$ ($\pm$1.96 SE))'))
      
      ylab_label_p2 <- ifelse(k == 1,
                              TeX(r'(RE($\hat{\beta}_0$))'),
                              TeX(r'(RE($\hat{\beta}_1$))'))
      
      p1 <- result %>%
        # filter(!(method %in% c("SL", "RF"))) %>% 
        ggplot(aes(x = method, y = est, color = method)) + 
        geom_errorbar(aes(ymin = lci, ymax = rci), width = .2, size = 1.2) +
        geom_point(size = 4) +
        theme_bw() +
        scale_color_manual(guide = "none", values = mypalette) +
        xlab(NULL) +
        ylab(ylab_label_p1) +
        theme(axis.text.x = element_text(size = 11, color = "black", face = "bold"),
              axis.text.y = element_text(size = 11),
              axis.title = element_text(size = 12),
              plot.title = element_text(size = 12, hjust = 0.5))
      if (k == 1) {
        p1 <- p1 + ggtitle(plot_title)
      }
      
      p2 <- result %>%
        # filter(!(method %in% c("SL", "RF"))) %>% 
        ggplot(aes(x = method, y = re, group = 1)) +
        geom_point(aes(color = method), size = 4) +
        # geom_line() +
        geom_hline(yintercept = 1, linetype = 2) +
        coord_cartesian(ylim = c(0.8, 2.3)) +
        #scale_y_continuous(n.breaks = 3) +
        theme_bw() +
        scale_color_manual(guide = "none", values = mypalette) +
        xlab(NULL) +
        ylab(ylab_label_p2) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 11))
      
      plot_est_collected <- c(plot_est_collected, list(p1))
      plot_re_collected <- c(plot_re_collected, list(p2))
    }
  }
}


# pdf(file = "analysis_result_plot (eif, continuous, cf, time).pdf", width = 16, height = 6)
# plot_est_collected[[1]] + plot_re_collected[[1]] + plot_est_collected[[2]] + plot_re_collected[[2]] +
#   plot_est_collected[[3]] + plot_re_collected[[3]] +
#   plot_layout(nrow = 2, ncol = 3, byrow = FALSE, height = rep(c(2, 1), 3))
# dev.off()

pdf(file = "analysis_result_plot (main paper).pdf", width = 16, height = 12)
plot_est_collected[[1]] + plot_est_collected[[4]] + plot_est_collected[[7]] +
  plot_re_collected[[1]] + plot_re_collected[[4]] + plot_re_collected[[7]] +
  plot_est_collected[[2]] + plot_est_collected[[5]] + plot_est_collected[[8]] + 
  plot_re_collected[[2]] + plot_re_collected[[5]] + plot_re_collected[[8]] +
  plot_est_collected[[3]] + plot_est_collected[[6]] + plot_est_collected[[9]] +
  plot_re_collected[[3]] + plot_re_collected[[6]] + plot_re_collected[[9]] +
  plot_layout(nrow = 6, ncol = 3, byrow = TRUE, height = rep(c(10, 5), 6))
dev.off()
