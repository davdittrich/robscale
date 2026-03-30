library(boot)
library(dplyr)
library(tidyr)
library(purrr)

#' Compute a memory- and CPU-aware worker count for mclapply.
#'
#' Reads /proc/meminfo (Linux) to estimate available RAM and current swap
#' pressure. Caps workers at the minimum of:
#'   1. BENCH_MAX_WORKERS env var (default 4) — hard ceiling
#'   2. floor(avail_ram * 0.80 / worker_mem_mb) — memory-based ceiling
#'   3. physical_cores - 1 — leave one core for the orchestration layer
#'
#' Emits INFO/WARNING messages so over-commitment is visible in the log.
#'
#' @param context   Label for log messages (e.g. "bootstrap CI")
#' @param worker_mem_mb  Estimated peak RSS per worker in MB (default 1024)
#' @return Integer worker count, at least 1
safe_worker_count <- function(context = "mclapply", worker_mem_mb = 1024L) {
  n_phys    <- parallel::detectCores(logical = FALSE)
  # Default: leave one core for the targets/callr orchestration layer.
  # Override with BENCH_MAX_WORKERS env var (e.g. for CI or restricted hosts).
  hard_cap  <- as.integer(Sys.getenv("BENCH_MAX_WORKERS", as.character(n_phys - 1L)))

  avail_mb  <- NA_real_
  swap_used_gb <- NA_real_

  if (file.exists("/proc/meminfo")) {
    lines    <- readLines("/proc/meminfo", warn = FALSE)
    parse_kb <- function(pat) {
      hit <- grep(pat, lines, value = TRUE)
      if (!length(hit)) return(NA_real_)
      as.numeric(gsub("[^0-9]", "", hit[1]))
    }
    avail_mb     <- parse_kb("^MemAvailable:") / 1024
    swap_total_kb <- parse_kb("^SwapTotal:")
    swap_free_kb  <- parse_kb("^SwapFree:")
    if (!is.na(swap_total_kb) && !is.na(swap_free_kb))
      swap_used_gb <- (swap_total_kb - swap_free_kb) / (1024^2)
  }

  # Memory-based ceiling: each worker needs ~worker_mem_mb; keep 5 % headroom.
  mem_cap <- if (!is.na(avail_mb)) {
    max(1L, floor(avail_mb * 0.95 / worker_mem_mb))
  } else {
    hard_cap
  }

  chosen <- max(1L, min(hard_cap, n_phys - 1L, mem_cap))

  # Always emit a line so resource decisions are visible in the pipeline log.
  message(sprintf(
    "INFO  [%s] workers=%d  (hard_cap=%d  mem_cap=%d  phys_cores=%d  avail_mb=%.0f)",
    context, chosen, hard_cap, mem_cap, n_phys,
    if (is.na(avail_mb)) -1 else avail_mb
  ))

  if (!is.na(swap_used_gb) && swap_used_gb > 1)
    message(sprintf(
      "WARN  [%s] %.1f GB swap already in use — results may be slower; workers reduced to %d.",
      context, swap_used_gb, chosen
    ))

  chosen
}

#' Compute bootstrap confidence interval for speedup
#'
#' @param target_times Vector of timings for the new implementation
#' @param ref_times    Vector of timings for the reference implementation
#' @param R            Number of bootstrap replicates
#' @param method       CI method: "percentile" (default) or "bca".
#'   Override via the BENCH_CI_METHOD environment variable.
compute_speedup_ci <- function(target_times, ref_times, R = 2000,
                               method = Sys.getenv("BENCH_CI_METHOD", "percentile")) {
  method <- match.arg(method, c("percentile", "bca"))

  len <- min(length(target_times), length(ref_times))
  boot_data <- data.frame(
    target = as.numeric(target_times[1:len]),
    ref    = as.numeric(ref_times[1:len])
  )

  ratio_median <- function(d, i) {
    median(d[i, "ref"]) / median(d[i, "target"])
  }

  # Bootstrap — run serially: this function is already called from mclapply,
  # so inner boot parallelism would nest forks (n_cores^2 processes, crash risk).
  set.seed(42)
  b <- boot::boot(boot_data, ratio_median, R = R, parallel = "no")

  # Degenerate fallback used by both paths when boot.ci returns NULL/short.
  degenerate_ci <- function() {
    if (length(b$t0) >= 1L && !is.na(b$t0[1L])) c(b$t0[1L], b$t0[1L])
    else c(NA_real_, NA_real_)
  }

  ci <- if (method == "percentile") {
    tryCatch(
      {
        res <- boot::boot.ci(b, type = "perc")
        pct <- res$percent
        if (is.null(pct) || length(pct) < 5L) stop("percentile returned NULL or short")
        c(pct[4L], pct[5L])
      },
      error = function(e) degenerate_ci()
    )
  } else {
    # BCa — try BCa first, fall back to percentile, then degenerate.
    # boot.ci warns (does not error) when all t are equal, returning NULL for
    # $bca — treat NULL as an error so the fallback chain fires.
    tryCatch(
      {
        res <- boot::boot.ci(b, type = "bca")
        bca <- res$bca
        if (is.null(bca) || length(bca) < 5L) stop("BCa returned NULL or short")
        c(bca[4L], bca[5L])
      },
      error = function(e) {
        tryCatch(
          {
            res <- boot::boot.ci(b, type = "perc")
            pct <- res$percent
            if (is.null(pct) || length(pct) < 5L) stop("percentile returned NULL")
            c(pct[4L], pct[5L])
          },
          error = function(e2) degenerate_ci()
        )
      }
    )
  }

  median_speedup <- median(boot_data$ref) / median(boot_data$target)
  ci_width_frac  <- (ci[2L] - ci[1L]) / median_speedup

  list(
    median_speedup = median_speedup,
    ci_low         = ci[1L],
    ci_high        = ci[2L],
    ci_method      = method,
    cv_target      = sd(boot_data$target) / mean(boot_data$target),
    cv_ref         = sd(boot_data$ref)    / mean(boot_data$ref),
    ci_width_frac  = ci_width_frac,
    reliable       = ci_width_frac < 0.5
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

    n_cores <- safe_worker_count("bootstrap CI")

    indices <- seq_len(nrow(df))
    ll <- parallel::mclapply(indices, function(i) {
      # Belt-and-suspenders: silence any inherited TBB threadpool inside forks.
      RcppParallel::setThreadOptions(numThreads = 1L)
      row <- df[i, ]
      target_data <- as.numeric(row[[target_col]][[1]])
      ref_data <- as.numeric(row[[ref_col]][[1]])

      analysis <- compute_speedup_ci(target_data, ref_data)
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
