#!/usr/bin/env Rscript
# Micro-benchmark: robScale and robLoc at small n to measure per-call overhead.
# Targets the stack frame penalty identified in the revss v3 comparison.
#
# Usage: Rscript benchmarks/small_n_stack_bench.R

library(bench)
library(robscale)

ns <- c(4L, 8L, 16L, 32L, 64L, 128L, 256L)

cat(sprintf("\n  Platform: %s (%s)\n", Sys.info()["sysname"], Sys.info()["machine"]))
cat(sprintf("  R: %s\n", R.version.string))
cat(sprintf("  robscale: %s\n\n", packageVersion("robscale")))

results <- vector("list", length(ns))
for (i in seq_along(ns)) {
  n <- ns[i]
  set.seed(42L)
  x <- rnorm(n)
  iters <- if (n <= 64) 10000L else 5000L
  bm <- bench::mark(
    robScale = robscale::robScale(x),
    robLoc   = robscale::robLoc(x),
    adm      = robscale::adm(x),
    mad_scaled = robscale::mad_scaled(x),
    check = FALSE, iterations = iters
  )
  bm <- as.data.frame(bm)[, c("expression", "median")]
  bm$expr <- as.character(bm$expression)
  bm$n <- n
  results[[i]] <- bm[, c("expr", "n", "median")]
}

all_bm <- do.call(rbind, results)

fmt <- function(t) {
  us <- as.numeric(t) * 1e6
  ifelse(us >= 1000, sprintf("%.1fms", us / 1000), sprintf("%.1fus", us))
}

cat("  n      robScale   robLoc     adm        mad_scaled\n")
cat("  ", strrep("-", 60), "\n")
for (n in ns) {
  sub <- all_bm[all_bm$n == n, ]
  rs  <- fmt(as.numeric(sub$median[sub$expr == "robScale"]))
  rl  <- fmt(as.numeric(sub$median[sub$expr == "robLoc"]))
  ad  <- fmt(as.numeric(sub$median[sub$expr == "adm"]))
  ms  <- fmt(as.numeric(sub$median[sub$expr == "mad_scaled"]))
  cat(sprintf("  %-6d %-10s %-10s %-10s %-10s\n", n, rs, rl, ad, ms))
}
cat("\n")
