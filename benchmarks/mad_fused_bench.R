#!/usr/bin/env Rscript
# P4 benchmark: fused median-then-MAD performance.
# Measures mad_scaled and ensemble at various n to detect memory-traffic improvement.

library(robscale)

ns <- c(4, 16, 64, 256, 1024, 4096, 16384, 65536)

get_iters <- function(n) {
  if (n <= 64) return(10000L)
  if (n <= 1024) return(5000L)
  if (n <= 16384) return(2000L)
  return(500L)
}

cat(sprintf("\n  Platform: %s (%s)\n", Sys.info()["sysname"], Sys.info()["machine"]))
cat(sprintf("  R: %s\n", R.version.string))
cat(sprintf("  robscale: %s\n\n", packageVersion("robscale")))

# Benchmark mad_scaled
cat(sprintf("  %-8s %-12s %-12s\n", "n", "mad_scaled", "robScale"))
cat("  ", strrep("-", 36), "\n", sep = "")

set.seed(42)
for (n in ns) {
  x <- rnorm(n)
  iters <- get_iters(n)

  # Warmup
  for (i in 1:min(100, iters)) { mad_scaled(x); robScale(x) }

  # mad_scaled timing
  t0 <- proc.time()["elapsed"]
  for (i in 1:iters) mad_scaled(x)
  t_mad <- (proc.time()["elapsed"] - t0) / iters * 1e6

  # robScale timing (also benefits from fused MAD in initial scale)
  t0 <- proc.time()["elapsed"]
  for (i in 1:iters) robScale(x)
  t_rs <- (proc.time()["elapsed"] - t0) / iters * 1e6

  cat(sprintf("  %-8d %-12s %-12s\n",
              n,
              sprintf("%.1fus", t_mad),
              sprintf("%.1fus", t_rs)))
}
