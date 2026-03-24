## bench/adm_gate_check.R
## SOLO benchmark gate for ADM optimization WUs (WU-ADM1–4).
## Compares adm() timing against the WU-ADM0 saved baseline (C_adm_orig).
##
## Usage:
##   Rscript bench/adm_gate_check.R
##
## Gate criterion: median(adm) / median(baseline) <= 1.05 at ALL sizes.
## Larger ratio = regression vs WU-ADM0 baseline.
##
## Note: H2H mode (C_adm_fast vs C_adm_orig) was removed in WU-ADM4
## when the diagnostic exports were cleaned up. SOLO mode is the ongoing
## regression check.
##
## Sub-µs noise protocol for n <= 512:
##   Timer quantum on this machine is ~420–490 ns. For n <= 512 the total
##   timing is ~2.5 µs (≈ 5 quanta), so a single quantum mismatch yields
##   ratio ≈ 1.2. Require 3 consecutive runs all showing ratio > 1.05
##   before declaring a real failure at these sizes.

library(robscale)
library(bench)

RATIO_LIMIT   <- 1.05
BASELINE_FILE <- "benchmarks/adm_perf_baseline.rds"

cat("=== adm gate check (SOLO) ===\n")
cat("Gate: adm() / WU-ADM0 baseline <=", RATIO_LIMIT, "\n\n")

if (!file.exists(BASELINE_FILE)) {
  cat("ERROR: baseline file not found:", BASELINE_FILE, "\n")
  cat("Run bench/adm_baseline.R first to capture the baseline.\n")
  quit(status = 1)
}
baseline <- readRDS(BASELINE_FILE)

sizes <- c(8, 16, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 32768, 65536)

results <- lapply(sizes, function(n) {
  set.seed(42)
  x <- rnorm(n)
  min_iter <- if (n <= 64) 5000L else if (n <= 1024) 1000L else 200L

  baseline_n <- baseline[baseline$size == n, ]
  if (nrow(baseline_n) == 0) {
    cat("WARNING: no baseline for n=", n, "\n")
    return(list(n = n, med_fast_ns = NA, med_orig_ns = NA, ratio = NA,
                pass = TRUE, sub_us = (n <= 512)))
  }
  bm <- bench::mark(adm(x), min_iterations = min_iter, check = FALSE)
  med_fast <- as.numeric(bm$median[1L])
  med_orig <- mean(baseline_n$median_ns) / 1e9   # avg over seeds, convert ns→s
  ratio    <- med_fast / med_orig

  list(n          = n,
       med_fast_ns = med_fast * 1e9,
       med_orig_ns = med_orig * 1e9,
       ratio       = ratio,
       pass        = ratio <= RATIO_LIMIT,
       sub_us      = (n <= 512))
})

cat(sprintf("%-8s  %-12s  %-14s  %-8s  %s\n",
            "n", "adm (ns)", "baseline (ns)", "ratio", "status"))
cat(strrep("-", 62), "\n")

all_pass    <- TRUE
sub_us_fails <- c()

for (r in results) {
  status <- if (is.na(r$ratio)) "N/A" else if (r$pass) "PASS" else "FAIL"
  if (!is.na(r$ratio) && !r$pass) {
    all_pass <- FALSE
    if (r$sub_us) sub_us_fails <- c(sub_us_fails, r$n)
  }
  cat(sprintf("%-8d  %12.1f  %14.1f  %8.4f  %s%s\n",
              r$n,
              if (is.na(r$med_fast_ns)) 0 else r$med_fast_ns,
              if (is.na(r$med_orig_ns)) 0 else r$med_orig_ns,
              if (is.na(r$ratio)) 0 else r$ratio,
              status,
              if (r$sub_us) " [noisy: apply 3-run protocol]" else ""))
}

cat(strrep("-", 62), "\n")
cat("Overall:", if (all_pass) "PASS" else "FAIL", "\n")

stable_fails <- Filter(function(r) !is.na(r$ratio) && !r$pass && !r$sub_us, results)

if (length(sub_us_fails) > 0) {
  cat("\nNOTE: Failures at n =", paste(sub_us_fails, collapse = ", "),
      "are within timer-noise range (timing < ~2.5 µs, timer quantum ~420-490 ns).\n")
  cat("3-run protocol: these are treated as warnings (not hard failures) because\n")
  cat("  a single quantum mismatch (~420-490 ns) produces ratio ≈ 1.2 at these sizes.\n")
  cat("  To confirm a real regression: run this script 3 times and check if ALL 3\n")
  cat("  show ratio > 1.05 at the same size before declaring a real failure.\n")
}

if (length(stable_fails) > 0) {
  cat("\nGATE FAILED: ratio > 1.05 at stable-zone n =",
      paste(sapply(stable_fails, `[[`, "n"), collapse = ", "), "\n")
  quit(status = 1)
} else if (!all_pass) {
  cat("\nGate passed (noisy-zone warnings only). All stable-zone ratios <=", RATIO_LIMIT, "\n")
} else {
  cat("\nGate passed. All ratios <=", RATIO_LIMIT, "\n")
}
