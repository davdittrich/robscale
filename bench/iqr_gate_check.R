# bench/iqr_gate_check.R
# Run after each Phase 1-7 code change.
# Compares iqr_scaled() performance against the Phase 0 baseline.
# Must use devtools::load_all() — not library(robscale).
#
# IMPORTANT — CPU frequency noise at n<=16:
#   At n=16 (~2µs), powersave governor causes bimodal quantization (~350ns spread)
#   that exceeds the 5% gate threshold (~97ns). Run via the wrapper to suppress it:
#     sudo bash bench/run_gate.sh
#   The wrapper sets governor=performance + FIFO-99 scheduling before running this script.
#
# Gate mode (auto-selected):
#   HEAD-TO-HEAD: if iqr_impl_orig() is exported in the current build, compare
#     iqr_impl() vs iqr_impl_orig() back-to-back in the SAME bench::mark() call.
#     This eliminates inter-session OS scheduling noise (bimodal ~3µs at n=16-17).
#   SAVED-BASELINE: otherwise, compare against benchmarks/iqr_perf_baseline.rds.
#
# Gate threshold (FIXED — never change):
#   ratio <= 1.05  for ALL sizes
# Noisy small-n measurements must be handled by better methodology (head-to-head +
# performance governor), not wider thresholds.
#
# Ratio aggregation: median across seeds per size (robust to per-seed variation).
library(bench)
devtools::load_all(quiet = TRUE)

# Check CPU governor — warn if not 'performance' (causes false failures at n<=16).
gov_file <- "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
if (file.exists(gov_file)) {
  gov <- readLines(gov_file, 1L)
  if (!identical(gov, "performance")) {
    cat(sprintf(paste0(
      "WARNING: CPU governor is '%s' (not 'performance').\n",
      "  Powersave frequency scaling causes ~350ns bimodal jitter at n<=16,\n",
      "  which exceeds the 5%% gate threshold (~97ns) and causes false failures.\n",
      "  For reliable n<=16 results, run via:\n",
      "    sudo bash bench/run_gate.sh\n\n"
    ), gov))
  }
}

K_IQR <- 0.741301109252801
sizes <- c(16L, 17L, 64L, 100L, 128L, 129L, 1000L, 2049L)
seeds <- c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L)
iters_for <- function(sz) if (sz <= 256L) 2000L else 500L

h2h_mode <- exists("iqr_impl_orig", mode = "function")
if (h2h_mode) {
  cat("Mode: HEAD-TO-HEAD (iqr_impl_orig detected — same-session comparison)\n\n")
} else {
  cat("Mode: SAVED-BASELINE (no iqr_impl_orig — comparing against saved RDS)\n\n")
  baseline <- readRDS("benchmarks/iqr_perf_baseline.rds")
  base_iqr <- baseline$iqr
  base_ens <- baseline$ens
}

# ---- measure ----
rows <- vector("list", length(sizes) * length(seeds))
idx  <- 1L
for (sz in sizes) {
  for (s in seeds) {
    set.seed(s)
    x <- rnorm(sz)
    if (h2h_mode) {
      bm <- bench::mark(
        orig = iqr_impl_orig(x, K_IQR),
        new  = iqr_impl(x, K_IQR),
        min_iterations = iters_for(sz),
        check = FALSE
      )
      expr_labels <- as.character(bm$expression)
      rows[[idx]] <- data.frame(
        size     = sz,
        seed     = s,
        base_ns  = as.numeric(bm$median[expr_labels == "orig"]) * 1e9,
        curr_ns  = as.numeric(bm$median[expr_labels == "new"])  * 1e9
      )
    } else {
      bm <- bench::mark(iqr_scaled(x), min_iterations = iters_for(sz), check = FALSE)
      br <- base_iqr[base_iqr$size == sz & base_iqr$seed == s, ]
      rows[[idx]] <- data.frame(
        size    = sz,
        seed    = s,
        base_ns = if (nrow(br) > 0) br$median_ns else NA_real_,
        curr_ns = as.numeric(bm$median) * 1e9
      )
    }
    idx <- idx + 1L
  }
}
results <- do.call(rbind, rows)

# ---- ensemble ----
set.seed(77)
x10 <- rnorm(10)
if (h2h_mode) {
  # Ensemble: no iqr_impl_orig shortcut; fall back to saved baseline for ensemble
  if (file.exists("benchmarks/iqr_perf_baseline.rds")) {
    base_ens_ns <- readRDS("benchmarks/iqr_perf_baseline.rds")$ens$median_ns
  } else {
    base_ens_ns <- NA_real_
  }
  bm_ens <- bench::mark(scale_robust(x10, n_boot = 50L),
                        min_iterations = 100L, check = FALSE)
  curr_ens_ns  <- as.numeric(bm_ens$median) * 1e9
} else {
  base_ens_ns <- base_ens$median_ns
  bm_ens      <- bench::mark(scale_robust(x10, n_boot = 50L),
                              min_iterations = 100L, check = FALSE)
  curr_ens_ns  <- as.numeric(bm_ens$median) * 1e9
}

# ---- per-seed detail ----
cat("=== Per-seed detail ===\n")
cat(sprintf("%-6s %-6s  %8s  %8s  %6s\n", "size", "seed", "base_ns", "curr_ns", "ratio"))
for (sz in sizes) {
  for (s in seeds) {
    r <- results[results$size == sz & results$seed == s, ]
    if (nrow(r) == 0 || is.na(r$base_ns)) next
    cat(sprintf("%-6d %-6d  %8.0f  %8.0f  %6.3f\n",
                sz, s, r$base_ns, r$curr_ns, r$curr_ns / r$base_ns))
  }
}

# ---- gate: median-of-seeds per size ----
cat("\n=== IQR performance gate (median-of-seeds) ===\n")
cat(sprintf("%-6s  %10s  %10s  %6s  %s\n",
            "size", "base_med_ns", "curr_med_ns", "ratio", "verdict"))

any_fail <- FALSE
for (sz in sizes) {
  thr    <- 1.05
  sub    <- results[results$size == sz, ]
  sub    <- sub[!is.na(sub$base_ns), ]
  base_m <- median(sub$base_ns)
  curr_m <- median(sub$curr_ns)
  ratio  <- curr_m / base_m
  verdict <- if (ratio <= thr) "PASS" else sprintf("FAIL (thr=%.2f)", thr)
  if (ratio > thr) any_fail <- TRUE
  cat(sprintf("%-6d  %10.0f  %10.0f  %6.3f  %s\n",
              sz, base_m, curr_m, ratio, verdict))
}

# Ensemble
if (!is.na(base_ens_ns)) {
  ens_ratio   <- curr_ens_ns / base_ens_ns
  ens_verdict <- if (ens_ratio <= 1.05) "PASS" else "FAIL (thr=1.05)"
  if (ens_ratio > 1.05) any_fail <- TRUE
  cat(sprintf("%-6s  %10.0f  %10.0f  %6.3f  %s  [ensemble n=10]\n",
              "ens", base_ens_ns, curr_ens_ns, ens_ratio, ens_verdict))
}

cat("\n")
if (any_fail) {
  cat("GATE FAILED — revert changes and investigate before proceeding.\n")
  quit(status = 1L)
} else {
  cat("GATE PASSED — all ratios within threshold.\n")
}
