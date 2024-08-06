rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

result <- readRDS("../simu-control_pattern/result_all.RDS")

library(tidyverse)


library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

# result$ML_method[is.na(result$ML_method)] <- ""
# result <- result %>% filter(generative_model != "complex")

const_S_halfrange <- 2

beta0_true <- 0.5
beta1_true <- 0.2

result <- result %>%
  mutate(beta0_in_ci = (beta0_uci >= beta0_true & beta0_lci <= beta0_true),
         beta1_in_ci = (beta1_uci >= beta1_true & beta1_lci <= beta1_true))

result_cp <- result %>%
  group_by(S_pattern,
           rand_prob_pattern,
           rand_prob_tuning_param,
           control_pattern,
           control_tuning_param,
           sample_size,
           total_T,
           moderator,
           control,
           estimator,
           ml_method,
           cross_fit
  ) %>%
  summarize(beta0_cp = mean(beta0_in_ci), beta1_cp = mean(beta1_in_ci),
            nsim_this_job = n()) %>%
  ungroup()

result_cp <- result_cp %>%
  mutate(beta0_cp_uci = beta0_cp + 1.96 * sqrt(beta0_cp * (1-beta0_cp) / nsim_this_job),
         beta0_cp_lci = beta0_cp - 1.96 * sqrt(beta0_cp * (1-beta0_cp) / nsim_this_job),
         beta1_cp_uci = beta1_cp + 1.96 * sqrt(beta1_cp * (1-beta1_cp) / nsim_this_job),
         beta1_cp_lci = beta1_cp - 1.96 * sqrt(beta1_cp * (1-beta1_cp) / nsim_this_job))

result_cp$cross_fit_label <- ifelse(result_cp$cross_fit, "cf", "")
result_cp$method_label <- factor(paste0(result_cp$estimator, ".", result_cp$ml_method,
                                        ".", result_cp$cross_fit_label))
result_cp$method_label <- fct_relevel(result_cp$method_label, 
                                      "eif_oracle..", 
                                      "wcls..",
                                      "eif.lm.",
                                      "eif.lm.cf",
                                      "eif.gam.",
                                      "eif.gam.cf",
                                      "eif.ranger.",
                                      "eif.ranger.cf",
                                      "eif.sl.smooth.",
                                      "eif.sl.smooth.cf",
                                      "eif.sl.all.",
                                      "eif.sl.all.cf")

methods_to_plot <- c("wcls..", "eif_oracle..", "eif.gam.", "eif.gam.cf",
                     "eif.ranger.", "eif.ranger.cf",
                     "eif.sl.all.", "eif.sl.all.cf")
result_cp <- result_cp %>% filter(sample_size == 100,
                                  control_tuning_param <= 3,
                                  method_label %in% methods_to_plot)

result_cp$method_plot <- result_cp$ml_method
result_cp$method_plot[result_cp$estimator == "eif_oracle"] <- "oracle"
result_cp$method_plot[result_cp$estimator == "wcls"] <- "wcls"
result_cp$method_plot <- fct_relevel(result_cp$method_plot, 
                                     "oracle",
                                     "gam",
                                     "ranger",
                                     "sl.all",
                                     "wcls")
result_cp$method_plot <- fct_recode(result_cp$method_plot,
                                    Oracle = "oracle",
                                    GAM = "gam",
                                    RF = "ranger",
                                    SL = "sl.all",
                                    WCLS = "wcls")

result_cp$linewidth_cf <- 1
result_cp$linewidth_cf[result_cp$cross_fit == TRUE] <- 1.2

result_cp$cross_fit <- ifelse(result_cp$cross_fit, "yes", "no")

result_cp <- result_cp %>% filter(sample_size == 100, control_tuning_param <= 1)

saveRDS(result_cp, "result_cp.RDS")


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

result <- readRDS("../simu-control_pattern/result_all.RDS")

library(tidyverse)


library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

# result$ML_method[is.na(result$ML_method)] <- ""
# result <- result %>% filter(generative_model != "complex")

const_S_halfrange <- 2

beta0_true <- 0.5
beta1_true <- 0.2

mse <- function(est, truth) {
  return(mean((est - truth)^2))
}

# select simulation - moderator = NULL

result <- result %>%
  filter(moderator == "NULL")

result_mse <- result %>%
  group_by(S_pattern,
           rand_prob_pattern,
           rand_prob_tuning_param,
           control_pattern,
           control_tuning_param,
           ar1_tuning_param,
           sample_size,
           total_T,
           moderator,
           control,
           estimator,
           ml_method,
           cross_fit
  ) %>%
  slice(1:100)

result_mse <- result_mse %>%
  summarize(beta0_mse = mse(beta0, beta0_true), beta1_mse = mse(beta1, beta1_true),
            nsim_this_job = n()) %>%
  ungroup()

result_mse$cross_fit_label <- ifelse(result_mse$cross_fit, "cf", "")
result_mse$method_label <- factor(paste0(result_mse$estimator, ".", result_mse$ml_method,
                                         ".", result_mse$cross_fit_label))
result_mse$method_label <- fct_relevel(result_mse$method_label, 
                                       "eif_oracle..", 
                                       "wcls..",
                                       "eif.lm.",
                                       "eif.lm.cf",
                                       "eif.gam.",
                                       "eif.gam.cf",
                                       "eif.ranger.",
                                       "eif.ranger.cf",
                                       "eif.sl.smooth.",
                                       "eif.sl.smooth.cf",
                                       "eif.sl.all.",
                                       "eif.sl.all.cf")

methods_to_plot <- c("wcls..", "eif_oracle..", "eif.gam.", "eif.gam.cf",
                     "eif.ranger.", "eif.ranger.cf",
                     "eif.sl.all.", "eif.sl.all.cf")
result_mse <- result_mse %>% filter(method_label %in% methods_to_plot)

result_mse$method_plot <- result_mse$ml_method
result_mse$method_plot[result_mse$estimator == "eif_oracle"] <- "oracle"
result_mse$method_plot[result_mse$estimator == "wcls"] <- "wcls"
result_mse$method_plot <- fct_relevel(result_mse$method_plot, 
                                      "oracle",
                                      "gam",
                                      "ranger",
                                      "sl.all",
                                      "wcls")
result_mse$method_plot <- fct_recode(result_mse$method_plot,
                                     Oracle = "oracle",
                                     GAM = "gam",
                                     RF = "ranger",
                                     SL = "sl.all",
                                     WCLS = "wcls")

result_mse$linewidth_cf <- 1
result_mse$linewidth_cf[result_mse$cross_fit == TRUE] <- 1.2

result_mse$cross_fit <- ifelse(result_mse$cross_fit, "yes", "no")

saveRDS(result_mse, "result_mse.RDS")

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

result <- readRDS("../simu-control_pattern/result_all.RDS")

library(tidyverse)


library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

# result$ML_method[is.na(result$ML_method)] <- ""
# result <- result %>% filter(generative_model != "complex")

const_S_halfrange <- 2

beta0_true <- 0.5
beta1_true <- 0.2

result_sd <- result %>%
  group_by(S_pattern,
           rand_prob_pattern,
           rand_prob_tuning_param,
           control_pattern,
           control_tuning_param,
           ar1_tuning_param,
           sample_size,
           total_T,
           moderator,
           control,
           estimator,
           ml_method,
           cross_fit
  ) %>%
  summarize(beta0_sd = sd(beta0), beta1_sd = sd(beta1),
            nsim_this_job = n()) %>%
  ungroup()

result_sd_nowcls <- result_sd %>% filter(estimator != "wcls")
result_sd_wcls <- result_sd %>% 
  filter(estimator == "wcls") %>%
  rename(beta0_sd_wcls = beta0_sd, beta1_sd_wcls = beta1_sd) %>%
  select(!c(estimator, ml_method, cross_fit))

result_re <- result_sd_nowcls %>%
  left_join(result_sd_wcls,
            by = c("S_pattern", "rand_prob_pattern", "rand_prob_tuning_param", 
                   "control_pattern", "control_tuning_param", 
                   "ar1_tuning_param",
                   "sample_size", "control", "moderator", "total_T")) %>%
  mutate(beta0_re = (beta0_sd_wcls / beta0_sd)^2, 
         beta1_re = (beta1_sd_wcls / beta1_sd)^2)

result_re$cross_fit_label <- ifelse(result_re$cross_fit, "cf", "")
result_re$method_label <- factor(paste0(result_re$estimator, ".", result_re$ml_method,
                                        ".", result_re$cross_fit_label))
result_re$method_label <- fct_relevel(result_re$method_label, 
                                      "eif_oracle..", 
                                      "wcls..",
                                      "eif.lm.",
                                      "eif.lm.cf",
                                      "eif.gam.",
                                      "eif.gam.cf",
                                      "eif.ranger.",
                                      "eif.ranger.cf",
                                      "eif.sl.smooth.",
                                      "eif.sl.smooth.cf",
                                      "eif.sl.all.",
                                      "eif.sl.all.cf")

methods_to_plot <- c("eif_oracle..", "eif.gam.", "eif.gam.cf",
                     "eif.ranger.", "eif.ranger.cf",
                     "eif.sl.all.", "eif.sl.all.cf")
result_re <- result_re %>% filter(sample_size == 100,
                                  method_label %in% methods_to_plot)

result_re$method_plot <- result_re$ml_method
result_re$method_plot[result_re$estimator == "eif_oracle"] <- "oracle"
result_re$method_plot <- fct_relevel(result_re$method_plot, 
                                     "oracle",
                                     "gam",
                                     "ranger",
                                     "sl.all")
result_re$method_plot <- fct_recode(result_re$method_plot,
                                    Oracle = "oracle",
                                    GAM = "gam",
                                    RF = "ranger",
                                    SL = "sl.all")


result_re$linewidth_cf <- 1
result_re$linewidth_cf[result_re$cross_fit == TRUE] <- 1.2

result_re$cross_fit <- ifelse(result_re$cross_fit, "yes", "no")

saveRDS(result_re, "result_re.RDS")


rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(patchwork)
library(latex2exp)

library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

result_re <- readRDS("result_re.RDS")
result_cp <- readRDS("result_cp.RDS")
result_mse <- readRDS("result_mse.RDS")

# result_re <- result_re %>% filter(control_pattern != "linear")
# result_cp <- result_cp %>% filter(control_pattern != "linear")
# result_mse <- result_mse %>% filter(control_pattern != "linear")


result_mse$control_label <- factor(result_mse$control)
levels(result_mse$control_label) <- c("dp" = TeX(r'(control = $t$)'),
                                      "dp S" = TeX(r'(control = $(t, Z_t)$)'))
result_mse$control_pattern_label <- factor(result_mse$control_pattern) %>%
    fct_recode("linear" = "linear",
               "simple nonlinear" = "dbeta",
               "periodic" = "sine",
               "step" = "step") %>%
    fct_relevel("linear", "simple nonlinear", "periodic", "step")

result_cp$control_label <- factor(result_cp$control)
levels(result_cp$control_label) <- c("dp" = TeX(r'(control = $t$)'),
                                     "dp S" = TeX(r'(control = $(t, Z_t)$)'))
result_cp$control_pattern_label <- factor(result_cp$control_pattern) %>%
    fct_recode("linear" = "linear",
               "simple nonlinear" = "dbeta",
               "periodic" = "sine",
               "step" = "step") %>%
    fct_relevel("linear", "simple nonlinear", "periodic", "step")

result_re$control_label <- factor(result_re$control)
levels(result_re$control_label) <- c("dp" = TeX(r'(control = $t$)'),
                                     "dp S" = TeX(r'(control = $(t, Z_t)$)'))
result_re$control_pattern_label <- factor(result_re$control_pattern) %>%
    fct_recode("linear" = "linear",
               "simple nonlinear" = "dbeta",
               "periodic" = "sine",
               "step" = "step") %>%
    fct_relevel("linear", "simple nonlinear", "periodic", "step")


mytheme <- theme(axis.text.x = element_text(size = 10),
                 axis.text.y = element_text(size = 10),
                 axis.title.x = element_text(size = 12),
                 axis.title.y = element_text(size = 12, face = "bold"),
                 strip.text = element_text(size = 13),
                 legend.title = element_text(size = 12, face = "bold"), 
                 legend.text = element_text(size = 11), 
                 plot.title = element_text(size = 14, hjust = 0.5),
                 panel.grid.minor = element_blank(),
                 strip.placement = "outside",
                 legend.margin = margin(t = -6),
                 strip.background = element_rect(color = "white", fill = "white")
                 )

plot_mse <- result_mse %>%
    filter(moderator == "NULL", control_tuning_param == 0.8) %>%
    ggplot(aes(x = sample_size, y = beta0_mse, 
               color = method_plot,
               linetype = cross_fit, linewidth = cross_fit)) +
    geom_line() +
    scale_linewidth_manual(name = "Cross-fitting", values = c(0.5, 1)) +
    facet_grid(control_pattern_label ~ control_label,
               labeller = labeller(control_label = label_parsed)) +
    scale_color_manual(name = "Estimator", values = mypalette) + 
    scale_linetype_discrete(name = "Cross-fitting") +
  coord_cartesian(ylim = c(0, 0.1)) +
    # scale_x_continuous(breaks = c(30, 50, 100)) +
    xlab("sample size") + 
    ylab("Mean Squared Error") +
    # ggtitle(TeX(r'(Mean Squared Error (Under $\lambda_1 = 1$))')) +
    theme_bw() + 
    mytheme +
    guides(color = guide_legend(override.aes = list(linewidth = 1)))

plot_cp <- result_cp %>%
    filter(moderator == "NULL") %>%
    ggplot(aes(x = control_tuning_param, y = beta0_cp, 
               color = method_plot,
               linetype = cross_fit, linewidth = cross_fit)) +
    geom_line() +
    scale_linewidth_manual(name = "Cross-fitting", values = c(0.5, 1)) +
    facet_grid(control_pattern_label ~ control_label,
               labeller = labeller(control_label = label_parsed)) +
    scale_color_manual(name = "Estimator", values = mypalette) + 
    scale_linetype_discrete(name = "Cross-fitting") +
    scale_y_continuous(breaks = c(0.85, 0.90, 0.95)) +
    xlab(TeX(r'($\lambda_1$ (nonlinearity in $\mu_t$))')) + 
    # ylab("CP") +
    ylab("Coverage Probability") + 
    theme_bw() +
    mytheme +
    guides(color = guide_legend(override.aes = list(linewidth = 1)))

plot_re <- result_re %>%
    filter(moderator == "NULL") %>%
    ggplot(aes(x = control_tuning_param, y = beta0_re, 
               color = method_plot,
               linetype = cross_fit, linewidth = cross_fit)) +
    # geom_jitter() + 
    geom_line() +
    scale_linewidth_manual(name = "Cross-fitting", values = c(0.5, 1)) +
    # geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    coord_cartesian(ylim = c(0.8, 1.8)) + 
    facet_grid(control_pattern_label ~ control_label,
               labeller = labeller(control_label = label_parsed)) +
    scale_color_manual(guide = "none", values = mypalette) + 
    scale_linetype_discrete(name = "Cross-fitting") +
    xlab(TeX(r'($\lambda_1$ (nonlinearity in $\mu_t$))')) + 
    # ylab("RE") +
    ylab("Relative Efficiency") +
    theme_bw() +
    mytheme

pdf(file = "simulation-eif-binary-control_pattern.pdf", width = 15/1.2, height = 8.5/1.2)
plot_mse + plot_cp + plot_re + 
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom')
dev.off()
