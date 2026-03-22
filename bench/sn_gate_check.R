# bench/sn_gate_check.R
# Run after each OPT-S1..S8 code change.
# Compares sn() performance against the pre-optimization baseline.
# Must use devtools::load_all() — not library(robscale).
#
# Gate mode (auto-selected):
#   HEAD-TO-HEAD (H2H): if C_sn_fast_orig() is exported in the current build,
#     compare C_sn_fast() vs C_sn_fast_orig() back-to-back in the SAME
#     bench::mark() call.  This eliminates inter-session OS scheduling noise.
#   SOLO: otherwise, capture baseline timings only.
#
# NOTE: Solo mode — no C_sn_fast_orig export. Cannot produce gate verdict.
# Captures baseline timings only.
#
# Gate threshold (FIXED — never change):
#   ratio <= 1.05  for ALL sizes
# Noisy small-n measurements must be handled by better methodology
# (head-to-head + performance governor), not wider thresholds.
#
# Sizes: 10, 16, 17, 50, 100, 128, 129, 500, 1000, 5000, 10000, 50000
library(bench)
tryCatch(devtools::load_all(quiet = TRUE), error = function(e) invisible(NULL))

# Check CPU governor — warn if not 'performance' (causes false failures at small n).
gov_file <- "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
if (file.exists(gov_file)) {
  gov <- readLines(gov_file, 1L)
  if (!identical(gov, "performance")) {
    cat(sprintf(paste0(
      "WARNING: CPU governor is '%s' (not 'performance').\n",
      "  Frequency scaling causes bimodal jitter at small n that can\n",
      "  exceed the 5%% gate threshold and cause false failures.\n",
      "  For reliable small-n results, run via:\n",
      "    sudo bash bench/run_gate.sh\n\n"
    ), gov))
  }
}

sizes <- c(10L, 16L, 17L, 50L, 100L, 128L, 129L, 500L, 1000L, 5000L, 10000L, 50000L)
seeds <- c(42L, 57L, 99L, 123L, 200L, 314L, 628L, 777L, 1024L, 1618L)

iters_for <- function(sz) {
  if (sz <= 17L)  100L
  else if (sz <= 500L) 100L
  else if (sz <= 5000L) 100L
  else 100L
}

h2h_mode <- exists("C_sn_fast_orig", mode = "function")

if (h2h_mode) {
  cat(paste0(
    "Mode: HEAD-TO-HEAD (C_sn_fast_orig detected)\n",
    "  Comparing C_sn_fast() vs C_sn_fast_orig() in same bench::mark() call.\n",
    "  Gate threshold: ratio <= 1.05 at ALL sizes.\n\n"
  ))
} else {
  cat(paste0(
    "Mode: SOLO (no C_sn_fast_orig export)\n",
    "# NOTE: Solo mode — no C_sn_fast_orig export. Cannot produce gate verdict.",
    " Captures baseline timings only.\n\n"
  ))
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
        orig = C_sn_fast_orig(x),
        new  = C_sn_fast(x),
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
      bm <- bench::mark(C_sn_fast(x), min_iterations = iters_for(sz), check = FALSE)
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

# ---- ensemble proxy ----
set.seed(77)
x_ens <- rnorm(10)
bm_ens <- bench::mark(scale_robust(x_ens, n_boot = 50L),
                      min_iterations = 100L, check = FALSE)
curr_ens_ns <- as.numeric(bm_ens$median) * 1e9

# ---- per-seed detail ----
cat("=== Per-seed detail ===\n")
if (h2h_mode) {
  cat(sprintf("%-7s %-6s  %8s  %8s  %6s\n", "size", "seed", "base_ns", "curr_ns", "ratio"))
} else {
  cat(sprintf("%-7s %-6s  %8s\n", "size", "seed", "curr_ns"))
}

for (sz in sizes) {
  for (s in seeds) {
    r <- results[results$size == sz & results$seed == s, ]
    if (nrow(r) == 0L) next
    if (h2h_mode) {
      if (!is.na(r$base_ns)) {
        cat(sprintf("%-7d %-6d  %8.0f  %8.0f  %6.3f\n",
                    sz, s, r$base_ns, r$curr_ns, r$curr_ns / r$base_ns))
      }
    } else {
      cat(sprintf("%-7d %-6d  %8.0f\n", sz, s, r$curr_ns))
    }
  }
}

# ---- gate / summary table ----
cat("\n=== Sn performance summary (median-of-seeds per size) ===\n")
thr <- 1.05

if (h2h_mode) {
  cat(sprintf("%-7s  %10s  %10s  %6s  %s\n",
              "size", "base_med_ns", "curr_med_ns", "ratio", "verdict"))
} else {
  cat(sprintf("%-7s  %10s\n", "size", "curr_med_ns"))
}

any_fail <- FALSE

for (sz in sizes) {
  sub    <- results[results$size == sz, ]
  curr_m <- median(sub$curr_ns)
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

# Ensemble row
if (h2h_mode) {
  # No orig ensemble to compare against in H2H mode — report solo
  cat(sprintf("%-7s  %10s  %10.0f  %6s  %s\n",
              "ens", "n/a", curr_ens_ns, "n/a", "[solo — no ens orig]"))
} else {
  cat(sprintf("%-7s  %10.0f\n", "ens", curr_ens_ns))
}

cat("\n")
if (h2h_mode) {
  if (any_fail) {
    cat("GATE FAILED — revert changes and investigate before proceeding.\n")
    quit(status = 1L)
  } else {
    cat("GATE PASSED — all ratios within threshold.\n")
  }
} else {
  cat("Solo mode — baseline timings captured. No gate verdict.\n")
  cat("Add C_sn_fast_orig export to activate H2H mode.\n")
}
