#!/usr/bin/env Rscript
# Phase 1: M-scale iteration count and Aitken fire-rate diagnostics
# Answers: at n=500-1000, how many rho_sum evaluations does robscale do?
# And how often does Aitken acceleration actually fire?
#
# Usage: Rscript bench/iteration_diagnostics.R

library(robscale)

n_vals  <- c(10, 20, 32, 64, 100, 200, 500, 1000, 2000, 5000)
n_reps  <- 200L  # random samples per n
set.seed(42)

cat(sprintf("%-6s  %7s  %6s  %6s  %8s  %8s  %10s\n",
    "n", "mean_oi", "sd_oi", "max_oi", "ak_rate", "rho_mean", "converge%"))
cat(strrep("-", 65), "\n")

results <- list()

for (n in n_vals) {
  outer_iters  <- integer(n_reps)
  aitken_fires <- integer(n_reps)
  rho_evals    <- integer(n_reps)
  converged    <- logical(n_reps)

  for (i in seq_len(n_reps)) {
    x <- rnorm(n)
    d <- robscale:::rob_scale_diag_impl(x)
    outer_iters[i]  <- d$outer_iters
    aitken_fires[i] <- d$aitken_fires
    rho_evals[i]    <- d$rho_evals
    converged[i]    <- d$converged
  }

  # Aitken rate = aitken_fires / outer_iters (avoid div-by-zero)
  ak_rate <- ifelse(outer_iters > 0, aitken_fires / outer_iters, 0)

  cat(sprintf("%-6d  %7.2f  %6.2f  %6d  %8.3f  %8.2f  %10.1f%%\n",
      n,
      mean(outer_iters),
      sd(outer_iters),
      max(outer_iters),
      mean(ak_rate),
      mean(rho_evals),
      100 * mean(converged)))

  results[[as.character(n)]] <- list(
    n            = n,
    outer_iters  = outer_iters,
    aitken_fires = aitken_fires,
    rho_evals    = rho_evals,
    converged    = converged,
    ak_rate      = ak_rate
  )
}

# Save for findings.md update
saveRDS(results, "benchmarks/iteration_diagnostics.rds")
cat("\nSaved to benchmarks/iteration_diagnostics.rds\n")

# Highlight the n=500-1000 regime
cat("\n=== Key finding for n=500-1000 regression analysis ===\n")
for (n in c(500, 1000)) {
  r <- results[[as.character(n)]]
  cat(sprintf("n=%d: mean rho_evals=%.1f  mean aitken_rate=%.1f%%  max_iters=%d\n",
      n, mean(r$rho_evals), 100*mean(r$ak_rate), max(r$outer_iters)))
}
