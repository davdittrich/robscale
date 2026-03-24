## bench/rob_scale_median_threshold.R
## WU-RS3: micro-benchmark the is_small median threshold for rob_scale_core.
##
## rob_scale_core calls median selection TWICE per invocation:
##   1. median(x)  → for centre estimate t
##   2. median(|x - t|) → for MAD/scale initialisation s_init
##
## We benchmark 2 consecutive calls of each path to mirror this pattern.
## Crossover criterion: n where median_net loses to fr_select (ratio > 1.0).
## Decision rule:
##   crossover in [11, 19]  → threshold confirmed at 16, no code change
##   crossover <= 10 or >= 20 → update rob_scale_core to separate threshold

library(robscale)
library(bench)

cat("=== WU-RS3: rob_scale_core is_small threshold calibration ===\n")
cat("ROBSCALE_SORT_MEDIAN_THRESHOLD = 16 (global constant)\n")
cat("Pattern: 2 consecutive calls (median + MAD)\n\n")

sizes <- c(8, 10, 12, 14, 15, 16, 17, 18, 20, 24, 32)
min_iter <- 5000L

results <- lapply(sizes, function(n) {
  set.seed(42); x <- rnorm(n)
  dev <- abs(x - median(x))  # pre-computed deviation vector for 2nd call

  bm <- bench::mark(
    net = {
      robscale:::bench_median_net_impl(x)
      robscale:::bench_median_net_impl(dev)
    },
    fr = {
      robscale:::bench_fr_select_impl(x)
      robscale:::bench_fr_select_impl(dev)
    },
    min_iterations = min_iter,
    max_time = 30,
    check = FALSE
  )

  med_net <- as.numeric(bm$median[1]) * 1e9  # ns
  med_fr  <- as.numeric(bm$median[2]) * 1e9
  ratio   <- med_net / med_fr  # > 1: net is slower than fr → fr wins → crossover

  list(n = n, med_net_ns = med_net, med_fr_ns = med_fr, ratio = ratio)
})

cat(sprintf("%-4s  %-12s  %-12s  %-8s  %s\n",
            "n", "2×net (ns)", "2×fr  (ns)", "net/fr", "winner"))
cat(strrep("-", 54), "\n")

crossover <- NA_integer_
for (r in results) {
  winner <- if (r$ratio <= 1.0) "net" else "fr"
  marker <- if (!is.na(crossover)) "" else if (r$ratio > 1.0) {
    crossover <<- r$n
    " ← crossover"
  } else ""
  cat(sprintf("%-4d  %12.1f  %12.1f  %8.4f  %s%s\n",
              r$n, r$med_net_ns, r$med_fr_ns, r$ratio, winner, marker))
}
cat(strrep("-", 54), "\n")

# Decision
cat("\n=== Decision ===\n")
if (is.na(crossover)) {
  cat("net wins at ALL tested sizes (crossover > 32)\n")
  cat("Result: threshold possibly too low — investigate larger n\n")
} else if (crossover >= 11 && crossover <= 19) {
  cat(sprintf("Crossover at n=%d — within [11, 19] of current threshold=16\n", crossover))
  cat("Result: THRESHOLD CONFIRMED ADEQUATE — no code change\n")
} else {
  cat(sprintf("Crossover at n=%d — OUTSIDE [11, 19]\n", crossover))
  if (crossover <= 10) {
    cat("Result: threshold should be LOWERED — fr wins earlier than expected\n")
  } else {
    cat("Result: threshold should be RAISED — net wins longer than expected\n")
  }
  cat("Action: introduce ROBSCALE_CORE_MEDIAN_THRESHOLD in robscale_config.h\n")
}
