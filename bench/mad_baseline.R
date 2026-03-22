# bench/mad_baseline.R
# Pre-optimization MAD performance baseline.
# Run ONCE, before Phase 2. Saved result used for regression gates in
# Phases 2, 3, 4a, 4b, 6, and 8. Phase 4b also saves a separate
# mad_from_data baseline used only by Phase 5.
# Must use devtools::load_all() — not library(robscale).
library(bench)
devtools::load_all(quiet = TRUE)

sizes <- c(64L, 65L, 100L, 1000L)
seeds <- c(42L, 99L, 123L)

rows <- vector("list", length(sizes) * length(seeds))
idx  <- 1L
for (sz in sizes) {
  for (s in seeds) {
    set.seed(s)
    x  <- rnorm(sz)
    bm <- bench::mark(mad_scaled(x), min_iterations = 500L, check = FALSE)
    rows[[idx]] <- data.frame(
      size      = sz,
      seed      = s,
      median_ns = as.numeric(bm$median) * 1e9
    )
    idx <- idx + 1L
  }
}
baseline_mad <- do.call(rbind, rows)

# Ensemble path: n=10 forces 7-estimator ensemble (below auto_switch threshold of 20)
set.seed(77)
x10     <- rnorm(10)
bm_ens  <- bench::mark(scale_robust(x10, n_boot = 50L),
                       min_iterations = 100L, check = FALSE)
baseline_ens <- data.frame(
  size      = 10L,
  n_boot    = 50L,
  median_ns = as.numeric(bm_ens$median) * 1e9
)

saveRDS(list(mad = baseline_mad, ens = baseline_ens),
        "benchmarks/mad_perf_baseline.rds",
        version = 2)

print(baseline_mad)
print(baseline_ens)
