# collect_results.R
#
# Aggregate result_tmp/<itask>.RDS into a tidy data frame in result_collected/.
# Sim 5.5 estimates only beta_0 (marginal CEE), so only beta_0 columns are kept.

suppressPackageStartupMessages(library(tidyverse))

# Rebuild the simulation design (must match simulation_R_code.R)

nsim_per_seed <- 100
n_seeds <- 10

setting_grid <- expand.grid(
    sample_size = c(30, 100),
    total_T     = c(10, 30, 50, 100, 200),
    stringsAsFactors = FALSE
) %>%
    mutate(setting_id = row_number())

simulation_design <- expand.grid(
    seed = 1:n_seeds,
    setting_id = setting_grid$setting_id,
    stringsAsFactors = FALSE
) %>%
    left_join(setting_grid, by = "setting_id") %>%
    arrange(setting_id, seed)

n_tasks <- nrow(simulation_design)
cat(sprintf("Total tasks in design: %d\n", n_tasks))

missing <- which(!file.exists(paste0("result_tmp/", 1:n_tasks, ".RDS")))
if (length(missing) > 0) {
    cat(sprintf("WARNING: %d result_tmp files are missing:\n", length(missing)))
    cat(paste0("  ", missing, "\n"))
    stop("Cannot collect results until all files are present.")
}
cat("All result_tmp files present.\n")

# Marginal CEE -> only beta_0 is estimated; eif_core / wcls_original return a
# 1-element beta_hat / beta_se / 1x2 conf_int.
extract_fit <- function(fit) {
    if (!is.null(fit$error_message)) {
        return(data.frame(beta0 = NA_real_,
                          beta0_se = NA_real_,
                          beta0_lci = NA_real_,
                          beta0_uci = NA_real_))
    }
    bh <- fit$beta_hat
    bs <- fit$beta_se_adjusted
    if (is.null(bs)) bs <- fit$beta_se
    ci <- fit$conf_int_adjusted
    if (is.null(ci)) ci <- fit$conf_int
    data.frame(
        beta0     = unname(bh[1]),
        beta0_se  = unname(bs[1]),
        beta0_lci = unname(ci[1, 1]),
        beta0_uci = unname(ci[1, 2])
    )
}

method_extractors <- list(
    wcls         = extract_fit,
    gam_cf       = extract_fit,
    sl_smooth_cf = extract_fit,
    oracle       = extract_fit
)

all_rows <- list()
for (itask in 1:n_tasks) {
    if (itask %% 10 == 0) cat(sprintf("  collecting task %d/%d\n", itask, n_tasks))
    setting <- simulation_design[itask, ]
    reps <- readRDS(paste0("result_tmp/", itask, ".RDS"))
    for (isim in seq_along(reps)) {
        rep <- reps[[isim]]
        for (method_name in names(method_extractors)) {
            row_df <- method_extractors[[method_name]](rep[[method_name]])
            row_df$method <- method_name
            row_df$itask <- itask
            row_df$isim <- isim
            row_df$seed <- setting$seed
            row_df$setting_id <- setting$setting_id
            row_df$sample_size <- setting$sample_size
            row_df$total_T <- setting$total_T
            all_rows[[length(all_rows) + 1]] <- row_df
        }
    }
}

result_long <- bind_rows(all_rows)

dir.create("result_collected", showWarnings = FALSE)
saveRDS(result_long, file = "result_collected/result_long.RDS")
write_csv(result_long, file = "result_collected/result_long.csv")

cat(sprintf("\nDone. %d rows collected (%d tasks x %d reps x %d methods).\n",
            nrow(result_long), n_tasks, nsim_per_seed, length(method_extractors)))
