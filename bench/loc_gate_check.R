# bench/loc_gate_check.R
# Run after each OPT-RL1..RL4 code change.
# Compares robLoc performance against an appropriate baseline.
#
# Gate mode (auto-selected):
#   HEAD-TO-HEAD (H2H): if rob_loc_fast_orig() is exported (WU-RL1 through WU-RL3),
#     compare rob_loc_impl() vs rob_loc_fast_orig() back-to-back in the SAME
#     bench::mark() call. Eliminates inter-session OS scheduling noise.
#   SOLO-GATE: if bench/loc_perf_baseline.rds exists and rob_loc_fast_orig is absent
#     (WU-RL4 onward), compare against the saved baseline.
#   SOLO-CAPTURE: if neither baseline nor _orig exist, capture timings only.
#
# Gate threshold (FIXED — never change):
#   ratio <= 1.05  for ALL sizes
#
# Sub-µs noise protocol: for n <= 16 (sub-2µs), a single ratio in [0.67, 1.43]
# is noise. Require 3 consecutive runs agreeing on direction before declaring
# failure. Handled by caller; this script reports raw ratios.
# See gotcha-sn-002 and gotcha-bench-001.
#
# Sizes: 4, 5, 6, 7, 8, 10, 16, 32, 64, 100, 200, 500, 1000
#   4..7:  straddle has_scale minobs=4 boundary and below AVX2 threshold
#   8:     at AVX2 fused-kernel threshold (RL1 target)
#   16,32: within micro_buf (MICRO_SIZE=64)
#   64,65: straddle micro_buf boundary
#   100+:  stack-buf path
#
# has_scale path: gated separately at sizes n in {4, 8, 16, 64, 100}
# Ensemble proxy: scale_robust(rnorm(n), n_boot = 50) at n in {50, 100, 200}.

library(bench)
tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))

# Check CPU governor
gov_file <- "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
if (file.exists(gov_file)) {
  gov <- readLines(gov_file, 1L)
  if (!identical(gov, "performance")) {
    cat(sprintf(paste0(
      "WARNING: CPU governor is '%s' (not 'performance').\n",
      "  Frequency scaling causes bimodal jitter at small n.\n",
      "  For reliable small-n results, run via:\n",
      "    sudo bash bench/run_gate.sh\n\n"
    ), gov))
  }
}

tol_nr <- sqrt(.Machine$double.eps)

baseline_path <- "bench/loc_perf_baseline.rds"
h2h_mode  <- exists("rob_loc_fast_orig", envir=asNamespace("robscale"), mode="function")
solo_gate <- !h2h_mode && file.exists(baseline_path)

if (h2h_mode) {
  cat(paste0(
    "Mode: HEAD-TO-HEAD (rob_loc_fast_orig detected)\n",
    "  Comparing rob_loc_impl() vs rob_loc_fast_orig() in same bench::mark() call.\n",
    "  Gate threshold: ratio <= 1.05 at ALL sizes.\n\n"
  ))
} else if (solo_gate) {
  cat(paste0(
    "Mode: SOLO-GATE (comparing against bench/loc_perf_baseline.rds)\n",
    "  Gate threshold: ratio <= 1.05 at ALL sizes.\n\n"
  ))
} else {
  cat(paste0(
    "Mode: SOLO-CAPTURE (no baseline, no _orig)\n",
    "  Capturing timings only — no gate verdict produced.\n",
    "  Run again after WU-RL3 removes _orig to save baseline.\n\n"
  ))
}

sizes <- c(4L, 5L, 6L, 7L, 8L, 10L, 16L, 32L, 64L, 100L, 200L, 500L, 1000L)
seeds <- c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L)

iters_for <- function(sz) {
  if (sz <= 16L)   200L
  else if (sz <= 200L) 100L
  else 50L
}

if (solo_gate) baseline <- readRDS(baseline_path)

# ---- measure rob_loc_impl (no-scale path) ----
rows <- vector("list", length(sizes) * length(seeds))
idx  <- 1L

for (sz in sizes) {
  for (s in seeds) {
    set.seed(s)
    x <- rnorm(sz)

    if (h2h_mode) {
      bm <- bench::mark(
        orig = robscale:::rob_loc_fast_orig(x, FALSE, 0.0, 80L, tol_nr),
        new  = robscale:::rob_loc_impl(x, FALSE, 0.0, 80L, tol_nr),
        min_iterations = iters_for(sz),
        check = FALSE
      )
      expr_labels <- as.character(bm$expression)
      rows[[idx]] <- data.frame(
        size    = sz,
        seed    = s,
        base_ns = as.numeric(bm$median[expr_labels == "orig"]) * 1e9,
        curr_ns = as.numeric(bm$median[expr_labels == "new"])  * 1e9
      )
    } else {
      bm <- bench::mark(
        robscale:::rob_loc_impl(x, FALSE, 0.0, 80L, tol_nr),
        min_iterations = iters_for(sz),
        check = FALSE
      )
      rows[[idx]] <- data.frame(
        size    = sz,
        seed    = s,
        base_ns = NA_real_,
        curr_ns = as.numeric(bm$median) * 1e9
      )
    }
    idx <- idx + 1L
  }
}
results <- do.call(rbind, rows)

# ---- measure has_scale path at selected sizes ----
hs_sizes  <- c(4L, 8L, 16L, 64L, 100L)
rows_hs   <- vector("list", length(hs_sizes) * length(seeds))
idx_hs    <- 1L
scale_val <- 1.5

for (sz in hs_sizes) {
  for (s in seeds) {
    set.seed(s)
    x <- rnorm(sz)

    if (h2h_mode) {
      bm <- bench::mark(
        orig = robscale:::rob_loc_fast_orig(x, TRUE, scale_val, 80L, tol_nr),
        new  = robscale:::rob_loc_impl(x, TRUE, scale_val, 80L, tol_nr),
        min_iterations = iters_for(sz),
        check = FALSE
      )
      expr_labels <- as.character(bm$expression)
      rows_hs[[idx_hs]] <- data.frame(
        size    = sz,
        seed    = s,
        base_ns = as.numeric(bm$median[expr_labels == "orig"]) * 1e9,
        curr_ns = as.numeric(bm$median[expr_labels == "new"])  * 1e9
      )
    } else {
      bm <- bench::mark(
        robscale:::rob_loc_impl(x, TRUE, scale_val, 80L, tol_nr),
        min_iterations = iters_for(sz),
        check = FALSE
      )
      rows_hs[[idx_hs]] <- data.frame(
        size    = sz,
        seed    = s,
        base_ns = NA_real_,
        curr_ns = as.numeric(bm$median) * 1e9
      )
    }
    idx_hs <- idx_hs + 1L
  }
}
results_hs <- do.call(rbind, rows_hs)

# ---- ensemble proxy ----
ens_sizes   <- c(50L, 100L, 200L)
ens_results <- vector("list", length(ens_sizes))
for (i in seq_along(ens_sizes)) {
  set.seed(77L + i)
  x_ens <- rnorm(ens_sizes[i])
  bm_ens <- bench::mark(scale_robust(x_ens, n_boot = 50L),
                        min_iterations = 20L, check = FALSE)
  ens_results[[i]] <- data.frame(
    ens_n   = ens_sizes[i],
    curr_ns = as.numeric(bm_ens$median) * 1e9
  )
}
ens_df <- do.call(rbind, ens_results)

# ---- report no-scale path ----
cat("=== rob_loc_impl performance (no-scale, median-of-seeds per size) ===\n")
thr <- 1.05

if (h2h_mode || solo_gate) {
  cat(sprintf("%-7s  %10s  %10s  %6s  %s\n",
              "size", "base_med_ns", "curr_med_ns", "ratio", "verdict"))
} else {
  cat(sprintf("%-7s  %10s\n", "size", "curr_med_ns"))
}

any_fail <- FALSE

for (sz in sizes) {
  sub    <- results[results$size == sz, ]
  curr_m <- median(sub$curr_ns, na.rm = TRUE)

  if (h2h_mode) {
    sub_ok <- sub[!is.na(sub$base_ns), ]
    base_m <- median(sub_ok$base_ns)
    curr_m <- median(sub_ok$curr_ns)
    ratio  <- curr_m / base_m
    verdict <- if (ratio <= thr) "PASS" else sprintf("FAIL (thr=%.2f)", thr)
    if (ratio > thr) any_fail <- TRUE
    cat(sprintf("%-7d  %10.0f  %10.0f  %6.3f  %s\n",
                sz, base_m, curr_m, ratio, verdict))
  } else if (solo_gate) {
    base_row <- baseline[baseline$size == sz, ]
    if (nrow(base_row) == 0L) {
      cat(sprintf("%-7d  %10s  %10.0f  %6s  %s\n", sz, "N/A", curr_m, "N/A", "no baseline"))
    } else {
      base_m <- base_row$curr_med_ns
      ratio  <- curr_m / base_m
      verdict <- if (ratio <= thr) "PASS" else sprintf("FAIL (thr=%.2f)", thr)
      if (ratio > thr) any_fail <- TRUE
      cat(sprintf("%-7d  %10.0f  %10.0f  %6.3f  %s\n",
                  sz, base_m, curr_m, ratio, verdict))
    }
  } else {
    cat(sprintf("%-7d  %10.0f\n", sz, curr_m))
  }
}

# ---- report has_scale path ----
cat(sprintf("\n=== rob_loc_impl has_scale=TRUE (scale=%.1f, median-of-seeds) ===\n", scale_val))
if (h2h_mode) {
  cat(sprintf("%-7s  %10s  %10s  %6s  %s\n",
              "size", "base_med_ns", "curr_med_ns", "ratio", "verdict"))
} else {
  cat(sprintf("%-7s  %10s\n", "size", "curr_med_ns"))
}

for (sz in hs_sizes) {
  sub    <- results_hs[results_hs$size == sz, ]
  curr_m <- median(sub$curr_ns, na.rm = TRUE)

  if (h2h_mode) {
    sub_ok <- sub[!is.na(sub$base_ns), ]
    base_m <- median(sub_ok$base_ns)
    curr_m <- median(sub_ok$curr_ns)
    ratio  <- curr_m / base_m
    verdict <- if (ratio <= thr) "PASS" else sprintf("FAIL (thr=%.2f)", thr)
    if (ratio > thr) any_fail <- TRUE
    cat(sprintf("%-7d  %10.0f  %10.0f  %6.3f  %s\n",
                sz, base_m, curr_m, ratio, verdict))
  } else {
    cat(sprintf("%-7d  %10.0f\n", sz, curr_m))
  }
}

# ---- report ensemble proxy ----
cat("\n=== Ensemble proxy (scale_robust, n_boot=50) ===\n")
cat(sprintf("%-7s  %10s\n", "ens_n", "curr_med_ns"))
for (i in seq_len(nrow(ens_df))) {
  cat(sprintf("%-7d  %10.0f\n", ens_df$ens_n[i], ens_df$curr_ns[i]))
}

# ---- save baseline / gate verdict ----
cat("\n")
if (h2h_mode) {
  if (any_fail) {
    cat("GATE FAILED — revert changes and investigate before proceeding.\n")
    quit(status = 1L)
  } else {
    cat("GATE PASSED — all ratios within threshold.\n")
  }
} else if (solo_gate) {
  if (any_fail) {
    cat("GATE FAILED — revert changes and investigate before proceeding.\n")
    quit(status = 1L)
  } else {
    cat("GATE PASSED — all ratios within threshold.\n")
  }
} else {
  # SOLO-CAPTURE: save baseline for future SOLO-GATE runs
  baseline_new <- data.frame(
    size        = sizes,
    curr_med_ns = vapply(sizes, function(sz) {
      median(results$curr_ns[results$size == sz], na.rm = TRUE)
    }, numeric(1L))
  )
  saveRDS(baseline_new, baseline_path)
  cat(sprintf("Baseline saved to %s\n", baseline_path))
  cat("Solo mode — no gate verdict. Re-run after code changes to compare against baseline.\n")
}
