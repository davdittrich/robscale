# benchmarks/rigorous_verify.R
library(callr)
library(bench)
library(boot)
library(dplyr)
library(tidyr)
library(purrr)

# Source analysis functions for BCa
source("benchmarks/analyze_results.R")

#' Build robscale with a specific flag
build_robscale <- function(fast_mode = 0) {
  lib <- tempfile(paste0("lib_", if(fast_mode) "fast" else "slow", "_"))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("Building robscale (FAST=%d) into %s...\n", fast_mode, lib))
  
  envs <- c(ROBSCALE_FAST = as.character(fast_mode), 
            R_LIBS_USER = .libPaths()[1],
            R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  
  res <- system2(R.home("bin/R"), 
                args = c("CMD", "INSTALL", ".", paste0("--library=", lib), "--no-multiarch", "--no-test-load"),
                env = paste0(names(envs), "=", envs),
                stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(res, "status"))) {
    stop("Building failed: ", paste(res, collapse = "\n"))
  }
  lib
}

#' Run a single estimator benchmark
run_benchmark <- function(lib_fast, lib_slow, estimator, n, iterations = NULL, min_time = 2.0) {
  cat(sprintf("Verifying %s (n=%d)...\n", estimator, n))
  
  measure <- function(lib, n, estimator, min_time, iters) {
    callr::r(function(lib, n, estimator, min_time, iters) {
      withr::with_libpaths(new = lib, action = "prefix", {
        library(robscale)
        library(bench)
        set.seed(42 + n)
        x <- rnorm(n)
        func <- getExportedValue("robscale", estimator)
        
        # Use a high number of iterations or time to get enough samples for BCa
        bm <- bench::mark(
          est = func(x),
          check = FALSE,
          min_iterations = if(!is.null(iters)) iters else 100,
          min_time = min_time
        )
        as.numeric(bm$time[[1]])
      })
    }, args = list(lib = lib, n = n, estimator = estimator, min_time = min_time, iters = iterations))
  }
  
  times_fast <- measure(lib_fast, n, estimator, min_time, iterations)
  times_slow <- measure(lib_slow, n, estimator, min_time, iterations)
  
  # Analyze with BCa (Lower Replicate count for speed, but user asked for 'sufficient')
  # We use 500 bootstrap replicates as a balance
  analysis <- compute_bca_speedup(times_fast, times_slow, R = 500)
  
  data.frame(
    estimator = estimator,
    n = n,
    median_slow = median(times_slow),
    median_fast = median(times_fast),
    speedup = analysis$median_speedup,
    ci_low = analysis$ci_low,
    ci_high = analysis$ci_high,
    status = if(analysis$ci_low >= 0.99) "✅ PASS" else "❌ FAIL"
  )
}

# Main Execution Script
# To be called by the agent to verify specific estimators or ranges

verify_all <- function(lib_fast, lib_slow) {
  # All sample sizes from the target suite
  n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)
  n_large <- c(32768, 65536, 131072, 262144, 524288, 1048576, 10000000)
  
  estimators_small <- c("robLoc", "robScale", "adm")
  estimators_large <- c("qn", "sn")
  
  results <- list()
  
  # 1. Micro-scale test for all (very important for regressions)
  for (est in c(estimators_small, estimators_large)) {
    for (n in c(3, 4, 16, 64)) {
      results[[length(results) + 1]] <- run_benchmark(lib_fast, lib_slow, est, n)
    }
  }
  
  # 2. Targeted verification for estimators at their specific crossover points
  # (To be expanded as needed)
  
  bind_rows(results)
}
