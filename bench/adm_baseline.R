# bench/adm_baseline.R
# Pre-optimization ADM performance baseline.
# Run ONCE before WU-ADM1. Saved result used for SOLO regression gates WU-ADM1–4.
# Must use devtools::load_all() — not library(robscale).
library(bench)
devtools::load_all(quiet = TRUE)

sizes  <- c(8L, 16L, 64L, 128L, 256L, 512L, 1024L, 2048L, 4096L, 8192L, 32768L, 65536L)
seeds  <- c(42L, 99L, 123L)
rows   <- vector("list", length(sizes) * length(seeds))
idx    <- 1L

for (sz in sizes) {
  for (s in seeds) {
    set.seed(s); x <- rnorm(sz)
    min_iter <- if (sz <= 64) 5000L else if (sz <= 1024) 1000L else 200L
    bm <- bench::mark(robscale:::C_adm_orig(x),
                      min_iterations = min_iter,
                      check = FALSE)
    rows[[idx]] <- data.frame(
      size      = sz,
      seed      = s,
      median_ns = as.numeric(bm$median[1L]) * 1e9
    )
    idx <- idx + 1L
  }
}

baseline_adm <- do.call(rbind, rows)
saveRDS(baseline_adm, "benchmarks/adm_perf_baseline.rds", version = 2)
print(baseline_adm)
cat("\nBaseline saved to benchmarks/adm_perf_baseline.rds\n")
