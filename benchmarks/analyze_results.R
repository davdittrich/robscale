library(boot)
library(dplyr)
library(tidyr)
library(purrr)

#' Compute BCa bootstrap confidence interval for speedup
#' @param target_times Vector of timings for the new implementation
#' @param ref_times Vector of timings for the reference implementation
#' @param R Number of bootstrap replicates
compute_bca_speedup <- function(target_times, ref_times, R = 500) {
  # Take min length to align samples
  len <- min(length(target_times), length(ref_times))
  boot_data <- data.frame(
    target = as.numeric(target_times[1:len]),
    ref = as.numeric(ref_times[1:len])
  )

  # Ratio of medians function
  ratio_median <- function(d, i) {
    target_med <- median(d[i, "target"])
    ref_med <- median(d[i, "ref"])
    ref_med / target_med
  }

  # Bootstrap
  set.seed(42)
  b <- boot::boot(boot_data, ratio_median, R = R, parallel = "no", ncpus = 1)

  # Try BCa
  ci <- tryCatch(
    {
      res <- boot::boot.ci(b, type = "bca")
      # bca is 4th and 5th element of $bca
      c(res$bca[4], res$bca[5])
    },
    error = function(e) {
      # Fallback to percentile if BCa fails (e.g. constant values)
      tryCatch(
        {
          res <- boot::boot.ci(b, type = "perc")
          c(res$percent[4], res$percent[5])
        },
        error = function(e2) {
          # If even percentile fails, check if all values are equal
          if (all(b$t == b$t0)) {
            return(c(b$t0, b$t0))
          }
          c(NA, NA)
        }
      )
    }
  )

  list(
    median_speedup = median(boot_data$ref) / median(boot_data$target),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

#' Analyze all benchmark results
analyze_benchmarks <- function(rob, legacy) {
  # Helper to flatten bench_mark objects
  flatten_bench <- function(bm, label) {
    if (is.null(bm)) {
      return(NULL)
    }
    bm %>%
      mutate(
        expr = as.character(expression),
        n = as.numeric(n),
        label = label
      ) %>%
      select(n, expr, time, label)
  }

  df_rob_m <- flatten_bench(rob$m_estimators, "robscale")
  df_rob_scale <- flatten_bench(rob$scale_estimators, "robscale")

  df_leg_revss <- flatten_bench(legacy$revss, "legacy")
  df_leg_robustbase <- flatten_bench(legacy$robustbase, "legacy")

  # 1. robscale vs revss (M-estimators)
  comp_leg_small <- df_rob_m %>%
    inner_join(df_leg_revss, by = c("n", "expr"), suffix = c("_rob", "_leg"))

  # 2. robscale vs robustbase (Scale estimators)
  comp_leg_large <- df_rob_scale %>%
    inner_join(df_leg_robustbase, by = c("n", "expr"), suffix = c("_rob", "_leg"))

  # Mapping function to apply bootstrap across rows
  run_analysis <- function(df, target_col, ref_col) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # Parallelize across rows using mclapply (base R parallel)
    # This is more efficient for our row-wise structure than inner boot parallelism
    n_cores <- parallel::detectCores(logical = FALSE)
    
    indices <- seq_len(nrow(df))
    ll <- parallel::mclapply(indices, function(i) {
      row <- df[i, ]
      target_data <- as.numeric(row[[target_col]][[1]])
      ref_data <- as.numeric(row[[ref_col]][[1]])
      
      analysis <- compute_bca_speedup(target_data, ref_data)
      cbind(row, as.data.frame(analysis))
    }, mc.cores = n_cores)
    
    bind_rows(ll)
  }

  # 3. New estimators: robscale vs existing R implementations
  # The expression names differ between robscale and legacy, so we join per-estimator
  df_rob_new <- flatten_bench(rob$new_estimators, "robscale")
  df_leg_new <- flatten_bench(legacy$new_estimators, "legacy")

  # Build comparison pairs: robscale vs each legacy competitor
  build_new_comp <- function(rob_df, leg_df, rob_expr, leg_expr, label) {
    r <- rob_df %>%
      filter(expr == rob_expr) %>%
      select(n, time_rob = time)
    l <- leg_df %>%
      filter(expr == leg_expr) %>%
      select(n, time_leg = time)
    inner_join(r, l, by = "n") %>% mutate(expr = label)
  }

  comp_new <- bind_rows(
    build_new_comp(df_rob_new, df_leg_new, "gmd", "gmd", "gmd vs Hmisc"),
    build_new_comp(df_rob_new, df_leg_new, "gmd", "gmd_gd", "gmd vs GiniDistance"),
    build_new_comp(df_rob_new, df_leg_new, "iqr_scaled", "iqr_scaled", "iqr_scaled vs stats"),
    build_new_comp(df_rob_new, df_leg_new, "iqr_scaled", "iqr_collapse", "iqr_scaled vs collapse"),
    build_new_comp(df_rob_new, df_leg_new, "mad_scaled", "mad_scaled", "mad_scaled vs stats"),
    build_new_comp(df_rob_new, df_leg_new, "mad_scaled", "mad_collapse", "mad_scaled vs collapse")
  )

  list(
    leg_small = run_analysis(comp_leg_small, "time_rob", "time_leg"),
    leg_large = run_analysis(comp_leg_large, "time_rob", "time_leg"),
    new_estimators = run_analysis(comp_new, "time_rob", "time_leg")
  )
}
