# bench/iqr_baseline.R
# Pre-optimization IQR performance baseline.
# Run ONCE, before Phase 1. Saved result is the reference for all regression
# gates in Phases 1 through 7.
# Must use devtools::load_all() — not library(robscale).
library(bench)
devtools::load_all(quiet = TRUE)

sizes <- c(16L, 17L, 64L, 100L, 128L, 129L, 1000L, 2049L)
seeds <- c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L)

# More iterations for sizes in the noisy sub-10µs range: 2000L for n≤256, 500L for n>256.
iters_for <- function(sz) if (sz <= 256L) 2000L else 500L

rows <- vector("list", length(sizes) * length(seeds))
idx  <- 1L
for (sz in sizes) {
  for (s in seeds) {
    set.seed(s)
    x  <- rnorm(sz)
    bm <- bench::mark(iqr_scaled(x), min_iterations = iters_for(sz), check = FALSE)
    rows[[idx]] <- data.frame(
      size      = sz,
      seed      = s,
      median_ns = as.numeric(bm$median) * 1e9
    )
    idx <- idx + 1L
  }
}
baseline_iqr <- do.call(rbind, rows)

# Ensemble proxy: n=10 forces 7-estimator ensemble (below auto_switch threshold of 20).
# IQR is one of 7 estimators called via compute_all_estimators() + ensemble_one_replicate().
set.seed(77)
x10     <- rnorm(10)
bm_ens  <- bench::mark(scale_robust(x10, n_boot = 50L),
                       min_iterations = 100L, check = FALSE)
baseline_ens <- data.frame(
  size      = 10L,
  n_boot    = 50L,
  median_ns = as.numeric(bm_ens$median) * 1e9
)

saveRDS(list(iqr = baseline_iqr, ens = baseline_ens),
        "benchmarks/iqr_perf_baseline.rds",
        version = 2)

print(baseline_iqr)
print(baseline_ens)
