#!/usr/bin/env Rscript
# Phase 2: Threshold crossover micro-benchmark
# Measures median_net vs std::nth_element (via floyd_rivest_select) for n=8..64.
# Determines empirical crossover n where FR beats median_net.
#
# Note: floyd_rivest_select for n < 600 calls std::nth_element internally.
# So this benchmark compares median_net (partial sort network) vs nth_element.
#
# Usage: Rscript bench/threshold_crossover_bench.R

library(robscale)
library(bench)

ns      <- c(8L, 10L, 12L, 14L, 16L, 18L, 20L, 24L, 28L, 32L, 40L, 48L, 56L, 64L)
# Loop-based timing: each function called many times inside one system.time()
# to overcome clock_gettime() overhead (~1 us) at tiny n.
n_iters <- 500000L
n_runs  <- 3L  # take minimum of 3 runs for stability

set.seed(42)

cat("Comparing median_net vs nth_element (via FR) for median extraction\n")
cat("Method: min of", n_runs, "runs x", n_iters, "iterations\n\n")

cat(sprintf("%-4s  %10s  %10s  %8s  %s\n",
    "n", "net_ns", "fr_ns", "ratio", "winner"))
cat(strrep("-", 50), "\n")

results <- list()

min_time <- function(fn, reps, n_runs) {
  times <- numeric(n_runs)
  for (r in seq_len(n_runs)) {
    times[r] <- system.time(for (i in seq_len(reps)) fn())["elapsed"]
  }
  min(times) / reps * 1e9  # ns per call
}

for (n in ns) {
  x <- rnorm(n)

  net_ns <- min_time(function() robscale:::bench_median_net_impl(x), n_iters, n_runs)
  fr_ns  <- min_time(function() robscale:::bench_fr_select_impl(x),  n_iters, n_runs)
  ratio  <- net_ns / fr_ns

  winner <- if (ratio > 1.05) "FR" else if (ratio < 0.95) "net" else "tie"

  cat(sprintf("%-4d  %10.1f  %10.1f  %8.3f  %s\n",
      n, net_ns, fr_ns, ratio, winner))

  results[[as.character(n)]] <- list(
    n       = n,
    net_ns  = net_ns,
    fr_ns   = fr_ns,
    ratio   = ratio,
    winner  = winner
  )
}

# Find crossover
winners <- sapply(results, `[[`, "winner")
crossover_candidates <- ns[winners == "FR"]
crossover_n <- if (length(crossover_candidates) > 0) min(crossover_candidates) else NA

cat("\n=== Summary ===\n")
cat("crossover_n =", crossover_n,
    "(first n where FR consistently beats median_net)\n")
cat("Recommendation for ROBSCALE_SORT_MEDIAN_THRESHOLD =",
    if (!is.na(crossover_n)) crossover_n - 2L else "needs manual inspection",
    "(one step below crossover for safety margin)\n")

# Save for findings.md update
saveRDS(list(results = results, crossover_n = crossover_n),
        "benchmarks/threshold_crossover.rds")
cat("\nSaved to benchmarks/threshold_crossover.rds\n")

# Also verify numerical agreement across threshold boundary
cat("\n=== Numerical agreement check (FR vs net, n=8..64) ===\n")
set.seed(99)
max_diff <- 0
for (n in ns) {
  for (rep in 1:20) {
    x <- rnorm(n)
    m1 <- robscale:::bench_median_net_impl(x)
    m2 <- robscale:::bench_fr_select_impl(x)
    d  <- abs(m1 - m2)
    if (d > max_diff) max_diff <- d
    if (d > 1e-13) {
      cat(sprintf("  DIFFER: n=%d rep=%d net=%.15g fr=%.15g diff=%.3e\n",
          n, rep, m1, m2, d))
    }
  }
}
cat(sprintf("Max absolute difference across all n=8..64: %.3e\n", max_diff))
cat("All within 1e-13:", max_diff < 1e-13, "\n")
