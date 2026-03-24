## bench/port2_gate.R
## WU-PORT-2: micro-buffer fix gate (mad_impl_auto, n=65..128 now uses stack path)
##
## CAPTURE mode: baseline measured with the fixed code (buf_micro[ROBSCALE_MICRO_BUFFER_SIZE])
## CHECK mode:   verify no regression vs baseline
##
## Stable zone (n=1000, n=5000): blocking — exit 1 on failure
## Noisy zone  (n=65, n=100, n=128): informational only

library(robscale)
library(bench)

RATIO_LIMIT   <- 1.05
BASELINE_FILE <- "bench/port2_pre.rds"
MODE <- if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
            toupper(commandArgs(trailingOnly = TRUE)[1]) == "CHECK") "CHECK" else "CAPTURE"

stable_sizes <- c(1000L, 5000L)
noisy_sizes  <- c(65L, 100L, 128L)

cat("=== port2 gate ===\n", "Mode:", MODE, "\n")

if (MODE == "CAPTURE") {
  all_sizes <- c(noisy_sizes, stable_sizes)
  pre <- lapply(all_sizes, function(n) {
    set.seed(42)
    x <- rnorm(n)
    bm <- bench::mark(robscale:::mad_impl_auto(x, 1), min_iterations = 5L, check = FALSE)
    med_ns <- as.numeric(bm$median[1]) * 1e9
    cat(sprintf("  n=%5d: %.0f ns\n", n, med_ns))
    list(n = n, median_ns = med_ns)
  })
  saveRDS(pre, BASELINE_FILE)
  cat("Baseline saved to", BASELINE_FILE, "\n")
} else {
  pre <- readRDS(BASELINE_FILE)
  stable_fails <- 0L
  for (b in pre) {
    set.seed(42)
    x <- rnorm(b$n)
    bm <- bench::mark(robscale:::mad_impl_auto(x, 1), min_iterations = 5L, check = FALSE)
    ratio <- as.numeric(bm$median[1]) * 1e9 / b$median_ns
    is_stable <- b$n > 512L
    pass <- ratio <= RATIO_LIMIT
    if (is_stable && !pass) stable_fails <- stable_fails + 1L
    cat(sprintf("n=%5d: %.4f %s %s\n", b$n, ratio,
                if (pass) "PASS" else "FAIL",
                if (is_stable) "(stable)" else "(noisy)"))
  }
  if (stable_fails > 0L) {
    cat("GATE FAILED\n")
    quit(status = 1)
  } else {
    cat("Gate passed.\n")
  }
}
