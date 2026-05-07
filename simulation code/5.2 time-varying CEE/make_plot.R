# make_plot.R
#
# Section 5.4 figure, mirroring the layout of Figure 1 in the existing
# manuscript: 3 plots (MSE, Coverage, Relative Efficiency) side-by-side, each
# with facet_grid(control_pattern ~ coef). Total 4 facet rows x 6 facet cols.
#
# Styling:
#   - RColorBrewer Set1 palette, fixed color per estimator across all panels
#   - Single bottom legend collected by patchwork
#   - RE panel suppresses its own legend (guide = "none")
#   - All methods solid line (we don't have a CF/non-CF distinction here)
#
# Usage:
#   Rscript make_plot.R

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(latex2exp)
    library(RColorBrewer)
})

# ---- True parameters (must match the dgm) ----
beta0_true <- 0.5
beta1_true <- -0.3

# ---- Load and add per-rep MSE / coverage indicators ----
result <- readRDS("result_collected/result_long.RDS")

result <- result %>%
    mutate(
        beta0_in_ci = (beta0_lci <= beta0_true) & (beta0_uci >= beta0_true),
        beta1_in_ci = (beta1_lci <= beta1_true) & (beta1_uci >= beta1_true),
        beta0_sq    = (beta0 - beta0_true)^2,
        beta1_sq    = (beta1 - beta1_true)^2
    )

# ---- Method factor with fixed colors (Set1; skip green slot from Fig 1) ----
method_levels <- c("oracle", "gam_cf", "sl_smooth_cf", "wcls")
method_labels <- c("Oracle", "GAM (CF)", "SL.smooth (CF)", "WCLS")

set1 <- brewer.pal(n = 5, "Set1")
mypalette <- setNames(set1[c(1, 2, 4, 5)], method_labels)

result <- result %>%
    mutate(method = factor(method, levels = method_levels, labels = method_labels))

# ---- Friendly control_pattern label, matching Figure 1 ----
pattern_recode <- c(linear = "linear",
                    dbeta  = "simple nonlinear",
                    sine   = "periodic",
                    step   = "step")
pattern_order <- c("linear", "simple nonlinear", "periodic", "step")

result <- result %>%
    mutate(control_pattern_label = factor(pattern_recode[control_pattern],
                                          levels = pattern_order))

# ---- Aggregate, then pivot to long on (beta_0, beta_1) ----
agg <- result %>%
    group_by(setting_id, sample_size, control_pattern, control_pattern_label,
             control_pattern_tuning_param, method) %>%
    summarise(beta0_mse = mean(beta0_sq, na.rm = TRUE),
              beta1_mse = mean(beta1_sq, na.rm = TRUE),
              beta0_cp  = mean(beta0_in_ci, na.rm = TRUE),
              beta1_cp  = mean(beta1_in_ci, na.rm = TRUE),
              .groups = "drop")

# RE = MSE_WCLS / MSE_method, joined within each setting
wcls_mse <- agg %>%
    filter(method == "WCLS") %>%
    select(setting_id,
           beta0_mse_wcls = beta0_mse,
           beta1_mse_wcls = beta1_mse)

agg <- agg %>%
    left_join(wcls_mse, by = "setting_id") %>%
    mutate(beta0_re = beta0_mse_wcls / beta0_mse,
           beta1_re = beta1_mse_wcls / beta1_mse)

# Pivot to long on coefficient
agg_long <- agg %>%
    pivot_longer(
        cols = matches("^beta[01]_(mse|cp|re)$"),
        names_to = c("coef", "metric"),
        names_pattern = "beta([01])_(mse|cp|re)",
        values_to = "value"
    ) %>%
    mutate(coef = factor(coef, levels = c("0", "1"),
                         labels = c("hat(beta)[0]", "hat(beta)[1]"))) %>%
    pivot_wider(names_from = metric, values_from = value)

# ---- Slices ----
agg_left   <- agg_long %>% filter(control_pattern_tuning_param == 1)
agg_middle <- agg_long %>% filter(sample_size == 100)

# ---- Theme matching Figure 1 ----
mytheme <- theme(axis.text.x      = element_text(size = 9),
                 axis.text.y      = element_text(size = 9),
                 axis.title.x     = element_text(size = 11),
                 axis.title.y     = element_text(size = 11, face = "bold"),
                 strip.text       = element_text(size = 11),
                 legend.title     = element_text(size = 11, face = "bold"),
                 legend.text      = element_text(size = 10),
                 plot.title       = element_text(size = 12, hjust = 0.5),
                 panel.grid.minor = element_blank(),
                 strip.placement  = "outside",
                 legend.margin    = margin(t = -6),
                 strip.background = element_rect(color = "white", fill = "white"))

facet_spec <- facet_grid(control_pattern_label ~ coef,
                         labeller = labeller(coef = label_parsed),
                         switch = "y")

# ---- The three plots ----
p_mse <- agg_left %>%
    ggplot(aes(x = sample_size, y = mse, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    facet_spec +
    scale_color_manual(name = "Estimator", values = mypalette) +
    xlab("sample size") +
    ylab("Mean Squared Error") +
    ggtitle(TeX(r'(MSE at $\lambda_1 = 1$)')) +
    theme_bw() + mytheme

p_cp <- agg_middle %>%
    ggplot(aes(x = control_pattern_tuning_param, y = cp, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    geom_hline(yintercept = 0.95, linetype = "dotted", alpha = 0.6) +
    facet_spec +
    scale_color_manual(name = "Estimator", values = mypalette) +
    coord_cartesian(ylim = c(0.85, 1)) +
    scale_y_continuous(breaks = c(0.85, 0.90, 0.95, 1.00)) +
    xlab(TeX(r'($\lambda_1$ (nonlinearity in $\mu_t$))')) +
    ylab("Coverage Probability") +
    ggtitle(TeX(r'(Coverage at $n = 100$)')) +
    theme_bw() + mytheme

p_re <- agg_middle %>%
    filter(method != "WCLS") %>%
    ggplot(aes(x = control_pattern_tuning_param, y = re, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.6) +
    facet_spec +
    scale_color_manual(values = mypalette, guide = "none") +
    xlab(TeX(r'($\lambda_1$ (nonlinearity in $\mu_t$))')) +
    ylab("Relative Efficiency") +
    ggtitle(TeX(r'(RE vs WCLS at $n = 100$)')) +
    theme_bw() + mytheme

# ---- Compose: 3 plots side-by-side, single bottom legend ----
combined <- (p_mse | p_cp | p_re) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

dir.create("plot", showWarnings = FALSE)
ggsave("plot/sim-5.4-time-varying-CEE.pdf", combined,
       width = 15 / 1.2, height = 8.5 / 1.2)

write_csv(agg, "result_collected/aggregated_summary.csv")

cat("Done. Figure: plot/sim-5.4-time-varying-CEE.pdf\n")
cat("Aggregated summary: result_collected/aggregated_summary.csv\n")
