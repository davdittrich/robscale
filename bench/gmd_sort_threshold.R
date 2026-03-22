#!/usr/bin/env Rscript
# OPT-G7: Benchmark full-sort crossover for ROBSCALE_SORT_MEDIAN_THRESHOLD
#
# Usage:
#   For each threshold T in {8, 10, 12, 14, 16}:
#     echo "PKG_CPPFLAGS = -DROBSCALE_SORT_MEDIAN_THRESHOLD=${T}" > src/Makevars.user
#     R CMD INSTALL . --no-test-load
#     Rscript bench/gmd_sort_threshold.R ${T}
#
# Decision rule: if any T < 16 reduces median time by >=10% for the majority
# of n <= 16, define ROBSCALE_GMD_SORT_THRESHOLD = T in robscale_config.h.
# Otherwise keep 16 and record "no improvement found".

args  <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1) args[1] else "default"

if (!requireNamespace("bench", quietly = TRUE))
  stop("Package 'bench' required: install.packages('bench')")

library(bench)
library(robscale)

ns <- c(4L, 6L, 8L, 10L, 12L, 14L, 16L, 20L, 32L)
cat(sprintf("# gmd_sort_threshold benchmark  threshold=%s  date=%s\n",
            label, Sys.Date()))
for (n in ns) {
  set.seed(1)
  x <- rnorm(n)
  bm <- bench::mark(gmd(x), iterations = 50000L, check = FALSE)
  cat(sprintf("threshold=%s n=%02d median=%.2f ns\n",
              label, n, as.numeric(bm$median) * 1e9))
}
