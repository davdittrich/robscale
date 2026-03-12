# benchmarks/run_micro_fast_vs_slow.R
library(callr)
library(bench)

# Function to build robscale with specific flag
build_robscale <- function(fast_mode = 0) {
  lib <- tempfile(paste0("lib_", if(fast_mode) "fast" else "slow", "_"))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  
  envs <- c(
    ROBSCALE_FAST = as.character(fast_mode),
    R_LIBS = .libPaths()[1],
    R_LIBS_USER = .libPaths()[1]
  )
  
  cat(sprintf("Building robscale (ROBSCALE_FAST=%d)...\n", fast_mode))
  res <- system2(R.home("bin/R"), 
                args = c("CMD", "INSTALL", ".", 
                        paste0("--library=", lib),
                        "--no-multiarch", "--no-test-load"),
                env = paste0(names(envs), "=", envs),
                stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(res, "status"))) {
    stop(sprintf("Installation failed (FAST=%d): ", fast_mode), paste(res, collapse = "\n"))
  }
  return(lib)
}

lib_slow <- build_robscale(0)
lib_fast <- build_robscale(1)

# Benchmark in isolated process
callr::r(function(lib_slow, lib_fast) {
  library(bench)
  
  measure <- function(lib, n, func_name) {
    callr::r(function(lib, n, func_name) {
      withr::with_libpaths(new = lib, action = "prefix", {
        library(robscale)
        set.seed(42 + n)
        x <- rnorm(n)
        func <- getExportedValue("robscale", func_name)
        res <- bench::mark(func(x), min_iterations = 100000, min_time = 1.0)
        return(as.numeric(res$median))
      })
    }, args = list(lib = lib, n = n, func_name = func_name))
  }
  
  compare <- function(n, func_name) {
    cat(sprintf("Benchmarking %s (n=%d)...\n", func_name, n))
    m_slow <- measure(lib_slow, n, func_name)
    m_fast <- measure(lib_fast, n, func_name)
    speedup <- m_slow / m_fast
    cat(sprintf("  Slow: %.2f µs\n  Fast: %.2f µs\n  Speedup: %.3fx\n", 
                m_slow * 1e6, m_fast * 1e6, speedup))
    if (speedup < 0.98) cat("  ❌ REGRESSION\n") else cat("  ✅ OK\n")
  }
  
  compare(4, "adm")
  compare(10, "adm")
  compare(32, "adm")
  compare(64, "adm")
  compare(5, "robLoc")
  compare(16, "robLoc")
  compare(64, "robLoc")
  compare(5, "robScale")
  compare(32, "robScale")
  compare(4, "qn")
  compare(16, "qn")
  compare(64, "qn")
  compare(4, "sn")
  compare(32, "sn")
  
}, args = list(lib_slow = lib_slow, lib_fast = lib_fast), show = TRUE)
