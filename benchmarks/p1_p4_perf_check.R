#!/usr/bin/env Rscript
# Performance check for P1 (incremental IQR), P2 (SIMD hints), P4 (fused MAD)
# Run against the INSTALLED version of robscale.
# Usage: Rscript --vanilla p1_p4_perf_check.R

library(robscale)
library(bench)

cat("=== System Info ===\n")
cat("Platform:", R.version$platform, "\n")
cat("R version:", R.version.string, "\n")
cat("robscale:", as.character(packageVersion("robscale")), "\n\n")

n_grid <- c(10, 16, 17, 32, 64, 128, 256, 512, 1000, 2048, 4096,
            10000, 50000, 100000, 500000, 1000000)

get_min_iters <- function(n) {
  if (n <= 128) 5000L
  else if (n <= 2048) 1000L
  else if (n <= 16384) 200L
  else if (n <= 100000) 50L
  else 10L
}

cat("=== IQR Benchmark (P1: incremental offset) ===\n")
iqr_results <- bench::press(
  n = n_grid,
  {
    set.seed(42 + n)
    x <- rnorm(n)
    bench::mark(
      robscale_iqr = robscale::iqr_scaled(x),
      base_iqr = IQR(x, type = 7) * 0.741301109252801,
      check = FALSE,
      min_iterations = get_min_iters(n),
      min_time = 0.5
    )
  }
)
cat("\nIQR results:\n")
iqr_summary <- iqr_results[, c("expression", "n", "median", "mem_alloc")]
print(iqr_summary, n = 100)

cat("\n=== MAD Benchmark (P4: fused median-then-MAD) ===\n")
mad_results <- bench::press(
  n = n_grid,
  {
    set.seed(42 + n)
    x <- rnorm(n)
    bench::mark(
      robscale_mad = robscale::mad_scaled(x),
      base_mad = stats::mad(x),
      check = FALSE,
      min_iterations = get_min_iters(n),
      min_time = 0.5
    )
  }
)
cat("\nMAD results:\n")
mad_summary <- mad_results[, c("expression", "n", "median", "mem_alloc")]
print(mad_summary, n = 100)

cat("\n=== GMD Benchmark (P2: SIMD reduction hints) ===\n")
gmd_results <- bench::press(
  n = n_grid,
  {
    set.seed(42 + n)
    x <- rnorm(n)
    bench::mark(
      robscale_gmd = robscale::gmd(x),
      check = FALSE,
      min_iterations = get_min_iters(n),
      min_time = 0.5
    )
  }
)
cat("\nGMD results:\n")
gmd_summary <- gmd_results[, c("expression", "n", "median", "mem_alloc")]
print(gmd_summary, n = 100)

cat("\n=== ADM Benchmark (P2: SIMD reduction hints) ===\n")
adm_results <- bench::press(
  n = n_grid,
  {
    set.seed(42 + n)
    x <- rnorm(n)
    bench::mark(
      robscale_adm = robscale::adm(x),
      check = FALSE,
      min_iterations = get_min_iters(n),
      min_time = 0.5
    )
  }
)
cat("\nADM results:\n")
adm_summary <- adm_results[, c("expression", "n", "median", "mem_alloc")]
print(adm_summary, n = 100)

cat("\n=== Ensemble Benchmark (all optimizations compound) ===\n")
ensemble_results <- bench::press(
  n = c(10, 32, 64, 128, 256, 512),
  {
    set.seed(42 + n)
    x <- rnorm(n)
    bench::mark(
      ensemble = robscale::scale_robust(x, method = "ensemble"),
      check = FALSE,
      min_iterations = 20L,
      min_time = 1.0
    )
  }
)
cat("\nEnsemble results:\n")
ens_summary <- ensemble_results[, c("expression", "n", "median", "mem_alloc")]
print(ens_summary, n = 100)

cat("\n=== IQR Memory Allocation Check ===\n")
for (nn in c(100, 10000, 1000000)) {
  set.seed(42)
  x <- rnorm(nn)
  m <- bench::mark(robscale::iqr_scaled(x), min_iterations = 10)
  cat(sprintf("n=%7d  mem_alloc=%s\n", nn, format(m$mem_alloc)))
}

cat("\n=== MAD Memory Allocation Check ===\n")
for (nn in c(100, 10000, 1000000)) {
  set.seed(42)
  x <- rnorm(nn)
  m <- bench::mark(robscale::mad_scaled(x), min_iterations = 10)
  cat(sprintf("n=%7d  mem_alloc=%s\n", nn, format(m$mem_alloc)))
}

cat("\nDone.\n")
