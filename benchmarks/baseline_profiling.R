#!/usr/bin/env Rscript
# Baseline Profiling — tracks robscale's absolute performance over time.
# Produces benchmarks/current_results.rds (list with $core and $ensemble).
#
# Reliability improvements (2026-03-21):
#   BM-1  callr isolation: runs in a fresh R subprocess to prevent library-path
#         contamination that previously caused 35x errors in stored baselines.
#   BM-2  3-seed pooling: time vectors from 3 independent random datasets are
#         concatenated before computing the summary statistics.  Reduces
#         sensitivity to one lucky/unlucky dataset without major runtime cost.
#   BM-3  Warmup + GC: 300-rep burn-in stabilises CPU P-states and instruction
#         cache; gc(full=TRUE) before each n group prevents GC pauses
#         contaminating the first iterations.
#   BM-7  min_time 0.5 → 1.0: more iterations per cell, narrower CI.
#
# Usage: Rscript benchmarks/baseline_profiling.R
# Expected runtime: ~10–20 min depending on hardware.

library(callr)
library(withr)

cat("Launching isolated measurement subprocess...\n")
cat(sprintf("Date: %s\n", Sys.Date()))

# Seeds are spaced 100 apart to avoid correlated .Random.seed states.
BENCH_SEEDS <- c(42L, 142L, 242L)

full_results <- callr::r(function(seeds) {
  library(robscale)
  library(bench)

  n_grid <- c(10L, 64L, 1024L, 16384L, 131072L)

  get_min_iters <- function(n) {
    if (n <= 64L)     10000L
    else if (n <= 1024L)  2000L
    else if (n <= 16384L)  500L
    else                    50L
  }

  # ── Helpers ───────────────────────────────────────────────────────────────

  # Pool time vectors from multiple bench::mark runs (one per seed).
  # Returns the first result with `time` list-column replaced by the pooled
  # concatenation across seeds, and summary stats recomputed from the pool.
  pool_bench_marks <- function(bm_list) {
    if (length(bm_list) == 1L) return(bm_list[[1L]])
    base <- bm_list[[1L]]
    for (i in seq_len(nrow(base))) {
      pooled_s <- unlist(lapply(bm_list, function(bm) as.numeric(bm$time[[i]])))
      base$time[[i]] <- pooled_s
    }
    # Recompute summary statistics from pooled times
    pool_s <- lapply(base$time, as.numeric)
    base$min         <- bench::as_bench_time(vapply(pool_s, min,    numeric(1L)))
    base$median      <- bench::as_bench_time(vapply(pool_s, median, numeric(1L)))
    base$mean        <- bench::as_bench_time(vapply(pool_s, mean,   numeric(1L)))
    base$max         <- bench::as_bench_time(vapply(pool_s, max,    numeric(1L)))
    base$n_itr       <- vapply(pool_s, length, integer(1L))
    base[["itr/sec"]]  <- 1 / vapply(pool_s, mean, numeric(1L))
    base$total_time  <- bench::as_bench_time(vapply(pool_s, sum,    numeric(1L)))
    base
  }

  # ── Warmup (BM-3) ─────────────────────────────────────────────────────────
  # Stabilise CPU P-states and L1/L2 cache before timed measurement.
  cat("  Warming up...\n")
  .wu_x <- rnorm(64L)
  for (.wu_i in seq_len(300L)) {
    robscale::mad_scaled(.wu_x)
    robscale::iqr_scaled(.wu_x)
    robscale::sn(.wu_x)
    robscale::qn(.wu_x)
    robscale::robScale(.wu_x)
    robscale::gmd(.wu_x)
    robscale::adm(.wu_x)
  }
  rm(.wu_x, .wu_i)
  gc(full = TRUE)

  # ── Core estimators (3-seed pooled) ───────────────────────────────────────
  results <- list()
  for (n in n_grid) {
    cat(sprintf("  Profiling n = %d ...\n", n))
    seed_bms <- lapply(seeds, function(seed) {
      set.seed(seed + n)
      x <- rnorm(n)
      gc(full = TRUE)         # flush GC before each configuration (BM-3)
      bench::mark(
        mad_scaled = robscale::mad_scaled(x),
        iqr_scaled = robscale::iqr_scaled(x),
        sn         = robscale::sn(x),
        qn         = robscale::qn(x),
        robScale   = robscale::robScale(x),
        gmd        = robscale::gmd(x),
        adm        = robscale::adm(x),
        check            = FALSE,
        min_iterations   = get_min_iters(n),
        min_time         = 1.0   # BM-7: was 0.5
      )
    })
    results[[as.character(n)]] <- pool_bench_marks(seed_bms)
  }

  # ── Ensemble (single-seed; too slow for 3-seed loop) ──────────────────────
  cat("  Profiling Ensemble...\n")
  ens_results <- list()
  for (n in c(10L, 64L, 1024L, 16384L)) {
    cat(sprintf("    n = %d (Ensemble)\n", n))
    set.seed(42L + n)
    x <- rnorm(n)
    gc(full = TRUE)
    bm <- bench::mark(
      ensemble = robscale::scale_robust(x, method = "ensemble"),
      check          = FALSE,
      min_iterations = 20L,
      min_time       = 1.0
    )
    bm$n <- n
    ens_results[[as.character(n)]] <- bm
  }

  list(
    core     = results,
    ensemble = ens_results,
    metadata = list(
      seeds          = seeds,
      n_grid         = n_grid,
      benchmark_date = as.character(Sys.Date()),
      r_version      = R.version$version.string,
      pkg_version    = as.character(packageVersion("robscale"))
    )
  )
}, args = list(seeds = BENCH_SEEDS), show = TRUE)

out_path <- "benchmarks/current_results.rds"
saveRDS(full_results, out_path)
cat(sprintf("\nBaseline results saved to %s\n", out_path))
cat(sprintf("  robscale %s | R %s | %d seeds | %s\n",
            full_results$metadata$pkg_version,
            R.version$major,
            length(BENCH_SEEDS),
            full_results$metadata$benchmark_date))
