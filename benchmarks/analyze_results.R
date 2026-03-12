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
  b <- boot::boot(boot_data, ratio_median, R = R)
  
  # Try BCa
  ci <- tryCatch({
    res <- boot::boot.ci(b, type = "bca")
    # bca is 4th and 5th element of $bca
    c(res$bca[4], res$bca[5])
  }, error = function(e) {
    # Fallback to percentile if BCa fails (e.g. constant values)
    tryCatch({
      res <- boot::boot.ci(b, type = "perc")
      c(res$percent[4], res$percent[5])
    }, error = function(e2) {
      # If even percentile fails, check if all values are equal
      if (all(b$t == b$t0)) {
        return(c(b$t0, b$t0))
      }
      c(NA, NA)
    })
  })
  
  list(
    median_speedup = median(boot_data$ref) / median(boot_data$target),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

#' Analyze all benchmark results
analyze_benchmarks <- function(rob_fast, rob_slow, legacy, rob_legacy_015 = NULL) {
  # Helper to flatten bench_mark objects
  flatten_bench <- function(bm, label) {
    if (is.null(bm)) return(NULL)
    bm %>%
      mutate(
        expr = as.character(expression),
        n = as.numeric(n),
        label = label
      ) %>%
      select(n, expr, time, label)
  }
  
  df_fast_m <- flatten_bench(rob_fast$m_estimators, "fast")
  df_fast_scale <- flatten_bench(rob_fast$scale_estimators, "fast")
  df_slow_m <- flatten_bench(rob_slow$m_estimators, "slow")
  df_slow_scale <- flatten_bench(rob_slow$scale_estimators, "slow")
  
  df_leg_revss <- flatten_bench(legacy$revss, "legacy")
  df_leg_robustbase <- flatten_bench(legacy$robustbase, "legacy")
  
  # 0.1.5 Gold Standard comparisons
  df_015_m <- if(!is.null(rob_legacy_015)) flatten_bench(rob_legacy_015$m_estimators, "0.1.5") else NULL
  df_015_scale <- if(!is.null(rob_legacy_015)) flatten_bench(rob_legacy_015$scale_estimators, "0.1.5") else NULL
  
  # Combine into comparison sets
  # 1. Optimized vs Legacy (External Packages)
  comp_leg_small <- df_fast_m %>%
    inner_join(df_leg_revss, by = c("n", "expr"), suffix = c("_fast", "_leg"))
  
  # 2. Optimized vs Legacy (External Packages - Large)
  comp_leg_large <- df_fast_scale %>%
    inner_join(df_leg_robustbase, by = c("n", "expr"), suffix = c("_fast", "_leg"))
    
  # 3. Optimized vs Unoptimized (Current Implementation)
  comp_fast_slow <- bind_rows(
    df_fast_m %>% inner_join(df_slow_m, by = c("n", "expr"), suffix = c("_fast", "_slow")),
    df_fast_scale %>% inner_join(df_slow_scale, by = c("n", "expr"), suffix = c("_fast", "_slow"))
  )

  # 4. Optimized vs 0.1.5 Gold Standard
  comp_fast_015 <- if(!is.null(df_015_m)) {
    bind_rows(
      df_fast_m %>% inner_join(df_015_m, by = c("n", "expr"), suffix = c("_fast", "_015")),
      df_fast_scale %>% inner_join(df_015_scale, by = c("n", "expr"), suffix = c("_fast", "_015"))
    )
  } else NULL
  
  # Mapping function to apply bootstrap across rows
  run_analysis <- function(df, target_col, ref_col) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      rowwise() %>%
      mutate(
        # In bench::mark, each entry in the 'time' list-column is a vector of timings.
        analysis = list(compute_bca_speedup(as.numeric(.data[[target_col]]), as.numeric(.data[[ref_col]])))
      ) %>%
      unnest_wider(analysis) %>%
      ungroup()
  }
  
  list(
    leg_small = run_analysis(comp_leg_small, "time_fast", "time_leg"),
    leg_large = run_analysis(comp_leg_large, "time_fast", "time_leg"),
    fast_slow = run_analysis(comp_fast_slow, "time_fast", "time_slow"),
    fast_015 = run_analysis(comp_fast_015, "time_fast", "time_015")
  )
}
