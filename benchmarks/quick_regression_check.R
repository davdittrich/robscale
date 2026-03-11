#!/usr/bin/env Rscript
# Quick regression check for robscale — lightweight standalone alternative to the
# full _targets.R pipeline. Runs three comparisons (FAST=1 vs legacy, FAST=0 vs
# legacy, FAST=1 vs FAST=0) on a reduced grid and prints pass/warn/fail verdicts.
#
# Usage: Rscript benchmarks/quick_regression_check.R          (all 3 comparisons, quick)
#        Rscript benchmarks/quick_regression_check.R --comp3  (FAST=1 vs FAST=0 only, full iterations)
# Expected runtime: ~3-5 min (default), ~5-10 min (--comp3)

library(bench)
library(withr)
library(dplyr, warn.conflicts = FALSE)

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
COMP3_ONLY <- "--comp3" %in% args

# ---------------------------------------------------------------------------
# Grid definitions
# ---------------------------------------------------------------------------
n_m <- c(4L, 8L, 16L, 64L, 256L, 1024L, 4096L, 16384L)
n_qnsn <- c(8L, 32L, 128L, 512L, 2048L, 8192L, 32768L, 131072L, 524288L, 1048576L)

# Quick mode: 10x fewer iterations than the full pipeline
get_iters_quick <- function(n) {
  if (n <= 128L)       200L
  else if (n <= 4096L) 100L
  else if (n <= 65536L) 50L
  else if (n <= 524288L) 30L
  else                  20L
}

# Full mode: matches the full _targets.R pipeline
get_iters_full <- function(n) {
  if (n <= 128L)        2000L
  else if (n <= 2048L)  1000L
  else if (n <= 16384L)  200L
  else if (n <= 1048576L) 50L
  else                    20L
}

get_iters <- if (COMP3_ONLY) get_iters_full else get_iters_quick

# ---------------------------------------------------------------------------
# Installation helpers (mirrors run_benchmarks.R)
# ---------------------------------------------------------------------------
install_robscale <- function(fast, lib_path) {
  dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)
  current_libs <- .libPaths()
  env_libs <- paste(current_libs, collapse = .Platform$path.sep)
  envs <- c(ROBSCALE_FAST = as.character(fast),
            R_LIBS = env_libs,
            R_LIBS_USER = env_libs)
  message("Installing robscale (FAST=", fast, ") into ", lib_path, " ...")
  res <- system2(
    R.home("bin/R"),
    args = c("CMD", "INSTALL", ".", paste0("--library=", lib_path),
             "--no-multiarch", "--no-test-load"),
    env = paste0(names(envs), "=", envs),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(res, "status"))) {
    stop("Installation failed (FAST=", fast, "): ", paste(res, collapse = "\n"))
  }
  invisible(res)
}

# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------
bench_m_estimators <- function(ns = n_m, loc_fn, scale_fn, adm_fn) {
  results <- vector("list", length(ns))
  for (i in seq_along(ns)) {
    n <- ns[i]
    set.seed(42L + n)
    x <- rnorm(n)
    iters <- get_iters(n)
    bm <- bench::mark(
      robLoc   = loc_fn(x),
      robScale = scale_fn(x),
      adm      = adm_fn(x),
      check = FALSE, iterations = iters
    )
    bm$n <- n
    results[[i]] <- bm
  }
  bind_rows(results)
}

bench_qnsn <- function(ns = n_qnsn, qn_call, sn_call) {
  results <- vector("list", length(ns))
  for (i in seq_along(ns)) {
    n <- ns[i]
    set.seed(42L + n)
    x <- rnorm(n)
    iters <- get_iters(n)
    bm <- bench::mark(
      qn = qn_call(x),
      sn = sn_call(x),
      check = FALSE, iterations = iters
    )
    bm$n <- n
    results[[i]] <- bm
  }
  bind_rows(results)
}

# ---------------------------------------------------------------------------
# Safely unload robscale if loaded
# ---------------------------------------------------------------------------
unload_robscale <- function() {
  if ("package:robscale" %in% search()) {
    detach("package:robscale", unload = TRUE, force = TRUE)
  }
}

# ---------------------------------------------------------------------------
# Run benchmarks for one robscale build
# ---------------------------------------------------------------------------
bench_robscale_build <- function(lib_path) {
  unload_robscale()
  withr::with_libpaths(new = lib_path, action = "prefix", {
    library(robscale)
    m <- bench_m_estimators(loc_fn   = robscale::robLoc,
                            scale_fn = robscale::robScale,
                            adm_fn   = robscale::adm)
    sc <- bench_qnsn(qn_call = robscale::qn,
                     sn_call = robscale::sn)
    list(m = m, sc = sc)
  })
}

# ---------------------------------------------------------------------------
# Run benchmarks for legacy packages
# ---------------------------------------------------------------------------
bench_legacy <- function() {
  library(revss)
  library(robustbase)
  m <- bench_m_estimators(loc_fn   = revss::robLoc,
                          scale_fn = revss::robScale,
                          adm_fn   = revss::adm)
  sc <- bench_qnsn(qn_call = robustbase::Qn,
                   sn_call = robustbase::Sn)
  list(m = m, sc = sc)
}

# ---------------------------------------------------------------------------
# Compute speedup table from two result sets
# ---------------------------------------------------------------------------
compute_speedup <- function(new_res, ref_res, label) {
  combine <- function(new_df, ref_df) {
    new_summ <- new_df |>
      mutate(expr = as.character(expression)) |>
      select(n, expr, median) |>
      rename(median_new = median)
    ref_summ <- ref_df |>
      mutate(expr = as.character(expression)) |>
      select(n, expr, median) |>
      rename(median_ref = median)
    inner_join(new_summ, ref_summ, by = c("n", "expr"))
  }
  merged <- bind_rows(
    combine(new_res$m, ref_res$m),
    combine(new_res$sc, ref_res$sc)
  ) |>
    mutate(
      speedup = as.numeric(median_ref) / as.numeric(median_new),
      status  = case_when(
        speedup < 0.95 ~ "FAIL",
        speedup < 1.00 ~ "WARN",
        TRUE           ~ "pass"
      ),
      comparison = label
    )
  merged
}

# ---------------------------------------------------------------------------
# Pretty-print a comparison table
# ---------------------------------------------------------------------------
print_table <- function(df, title) {
  cat("\n", strrep("=", 70), "\n")
  cat(" ", title, "\n")
  cat(strrep("=", 70), "\n")
  display <- df |>
    mutate(
      median_new_ms = sprintf("%.3f", as.numeric(median_new) * 1e3),
      median_ref_ms = sprintf("%.3f", as.numeric(median_ref) * 1e3),
      speedup_fmt   = sprintf("%.2fx", speedup)
    ) |>
    select(expr, n, median_new_ms, median_ref_ms, speedup_fmt, status)
  print(as.data.frame(display), right = FALSE, row.names = FALSE)
  n_fail <- sum(df$status == "FAIL")
  n_warn <- sum(df$status == "WARN")
  cat(sprintf("\n  Cells: %d pass, %d WARN, %d FAIL\n",
              nrow(df) - n_fail - n_warn, n_warn, n_fail))
}

# ===========================================================================
# Main
# ===========================================================================
main <- function() {
  start_time <- proc.time()

  if (COMP3_ONLY) {
    message("Mode: --comp3 (FAST=1 vs FAST=0 only, full iterations)")
  } else {
    message("Mode: default (all 3 comparisons, quick iterations)")
  }

  # 1. Install two builds
  lib_fast <- tempfile("robscale_fast_")
  lib_slow <- tempfile("robscale_slow_")
  install_robscale(fast = 1, lib_path = lib_fast)
  install_robscale(fast = 0, lib_path = lib_slow)

  # 2. Benchmark legacy packages (skip in --comp3 mode)
  if (!COMP3_ONLY) {
    message("\nBenchmarking legacy packages (revss + robustbase) ...")
    legacy <- bench_legacy()
  }

  # 3. Benchmark FAST=1
  message("\nBenchmarking robscale FAST=1 ...")
  fast <- bench_robscale_build(lib_fast)

  # 4. Benchmark FAST=0
  message("\nBenchmarking robscale FAST=0 ...")
  slow <- bench_robscale_build(lib_slow)

  # 5. Compute speedups
  if (!COMP3_ONLY) {
    comp1 <- compute_speedup(fast, legacy, "FAST1_vs_legacy")
    comp2 <- compute_speedup(slow, legacy, "FAST0_vs_legacy")
  }
  comp3 <- compute_speedup(fast, slow, "FAST1_vs_FAST0")

  # 6. Print results
  if (!COMP3_ONLY) {
    print_table(comp1, "Comparison 1: robscale FAST=1 vs legacy (revss/robustbase)")
    print_table(comp2, "Comparison 2: robscale FAST=0 vs legacy (revss/robustbase)")
  }
  print_table(comp3, "Comparison 3: robscale FAST=1 vs FAST=0")

  # 7. Overall verdict
  if (COMP3_ONLY) {
    all_results <- comp3
  } else {
    all_results <- bind_rows(comp1, comp2, comp3)
  }
  any_fail <- any(all_results$status == "FAIL")
  any_warn <- any(all_results$status == "WARN")

  cat("\n", strrep("=", 70), "\n")
  if (any_fail) {
    cat("  OVERALL VERDICT: FAIL\n")
  } else if (any_warn) {
    cat("  OVERALL VERDICT: PASS (with warnings)\n")
  } else {
    cat("  OVERALL VERDICT: PASS\n")
  }
  elapsed <- (proc.time() - start_time)[["elapsed"]]
  cat(sprintf("  Total time: %.1f seconds\n", elapsed))
  cat(strrep("=", 70), "\n")

  # 8. Save CSV
  suffix <- if (COMP3_ONLY) "_comp3" else ""
  out_file <- file.path("benchmarks",
                        paste0("quick_check_", format(Sys.Date(), "%Y%m%d"), suffix, ".csv"))
  write_df <- all_results |>
    mutate(
      median_new_s = as.numeric(median_new),
      median_ref_s = as.numeric(median_ref)
    ) |>
    select(comparison, expr, n, median_new_s, median_ref_s, speedup, status)
  write.csv(write_df, out_file, row.names = FALSE)
  message("Results saved to ", out_file)

  # Exit with non-zero status on failure
  if (any_fail) quit(status = 1)
}

main()
