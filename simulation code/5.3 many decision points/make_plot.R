# make_plot.R
#
# Section 5.5 figure: behavior under a large-T regime.
# Layout: 3 plots side-by-side (MSE, Coverage, RE-vs-WCLS), each with T on the
# x-axis. Lines colored by estimator; sample size n shown in facet rows.
# Mirrors the styling conventions of Figure 1 / sim 5.4.

suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(latex2exp)
    library(RColorBrewer)
})

beta0_true <- 0.5

result <- readRDS("result_collected/result_long.RDS")

result <- result %>%
    mutate(
        beta0_in_ci = (beta0_lci <= beta0_true) & (beta0_uci >= beta0_true),
        beta0_sq    = (beta0 - beta0_true)^2
    )

# Method factor with fixed colors (Set1; skip green slot for consistency
# with sim 5.4 / Figure 1 conventions)
method_levels <- c("oracle", "gam_cf", "sl_smooth_cf", "wcls")
method_labels <- c("Oracle", "GAM (CF)", "SL.smooth (CF)", "WCLS")

set1 <- brewer.pal(n = 5, "Set1")
mypalette <- setNames(set1[c(1, 2, 4, 5)], method_labels)

result <- result %>%
    mutate(method = factor(method, levels = method_levels, labels = method_labels),
           sample_size_label = factor(sample_size,
                                      levels = c(30, 100),
                                      labels = c("n == 30", "n == 100")))

agg <- result %>%
    group_by(setting_id, sample_size, sample_size_label, total_T, method) %>%
    summarise(beta0_mse = mean(beta0_sq, na.rm = TRUE),
              beta0_cp  = mean(beta0_in_ci, na.rm = TRUE),
              .groups = "drop")

wcls_mse <- agg %>%
    filter(method == "WCLS") %>%
    select(setting_id, beta0_mse_wcls = beta0_mse)

agg <- agg %>%
    left_join(wcls_mse, by = "setting_id") %>%
    mutate(beta0_re = beta0_mse_wcls / beta0_mse)

mytheme <- theme(axis.text.x      = element_text(size = 10),
                 axis.text.y      = element_text(size = 10),
                 axis.title.x     = element_text(size = 12),
                 axis.title.y     = element_text(size = 12, face = "bold"),
                 strip.text       = element_text(size = 12),
                 legend.title     = element_text(size = 12, face = "bold"),
                 legend.text      = element_text(size = 11),
                 plot.title       = element_text(size = 13, hjust = 0.5),
                 panel.grid.minor = element_blank(),
                 strip.placement  = "outside",
                 legend.margin    = margin(t = -6),
                 strip.background = element_rect(color = "white", fill = "white"))

facet_spec <- facet_wrap(~ sample_size_label, nrow = 2,
                         labeller = labeller(sample_size_label = label_parsed),
                         strip.position = "right")

p_mse <- agg %>%
    ggplot(aes(x = total_T, y = beta0_mse, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    facet_spec +
    scale_color_manual(name = "Estimator", values = mypalette) +
    scale_x_continuous(breaks = c(10, 30, 50, 100, 200)) +
    scale_y_log10() +
    xlab("T (number of decision points)") +
    ylab(TeX(r'(MSE($\hat{\beta}_0$))')) +
    ggtitle("Mean Squared Error") +
    theme_bw() + mytheme

p_cp <- agg %>%
    ggplot(aes(x = total_T, y = beta0_cp, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    geom_hline(yintercept = 0.95, linetype = "dotted", alpha = 0.6) +
    facet_spec +
    scale_color_manual(name = "Estimator", values = mypalette) +
    scale_x_continuous(breaks = c(10, 30, 50, 100, 200)) +
    coord_cartesian(ylim = c(0.85, 1)) +
    scale_y_continuous(breaks = c(0.85, 0.90, 0.95, 1.00)) +
    xlab("T (number of decision points)") +
    ylab("Coverage Probability") +
    ggtitle(TeX(r'(Coverage of $\hat{\beta}_0$)')) +
    theme_bw() + mytheme

p_re <- agg %>%
    filter(method != "WCLS") %>%
    ggplot(aes(x = total_T, y = beta0_re, color = method)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    geom_hline(yintercept = 1, linetype = "dotted", alpha = 0.6) +
    facet_spec +
    scale_color_manual(values = mypalette, guide = "none") +
    scale_x_continuous(breaks = c(10, 30, 50, 100, 200)) +
    xlab("T (number of decision points)") +
    ylab("Relative Efficiency") +
    ggtitle("RE vs WCLS") +
    theme_bw() + mytheme

combined <- (p_mse | p_cp | p_re) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

dir.create("plot", showWarnings = FALSE)
ggsave("plot/sim-5.5-large-T-regime.pdf", combined,
       width = 15 / 1.2, height = 6)

write_csv(agg, "result_collected/aggregated_summary.csv")

cat("Done. Figure: plot/sim-5.5-large-T-regime.pdf\n")
cat("Aggregated summary: result_collected/aggregated_summary.csv\n")
