#!/usr/bin/env Rscript
# adaptive_mad_ab.R — A/B: adaptive mad_scaled vs 0.2.0 (pure pdqselect) mad_scaled
#
# Loads each version from a separate library and benchmarks at sizes
# around the Zen 3 regression range.

suppressPackageStartupMessages(library(bench))

# --- load both versions via separate namespaces ---
baseline_env <- new.env(parent = baseenv())
adaptive_env <- new.env(parent = baseenv())

# Load baseline (pure pdqselect)
lib_b <- Sys.getenv("LIB_BASELINE", "~/lib_baseline")
lib_a <- Sys.getenv("LIB_ADAPTIVE", "~/lib_adaptive")

suppressPackageStartupMessages(
  library(robscale, lib.loc = lib_b, warn.conflicts = FALSE)
)
mad_baseline <- robscale::mad_scaled
cfg_baseline <- robscale:::get_qnsn_config()
unloadNamespace("robscale")

suppressPackageStartupMessages(
  library(robscale, lib.loc = lib_a, warn.conflicts = FALSE)
)
mad_adaptive <- robscale::mad_scaled
cfg_adaptive <- robscale:::get_qnsn_config()
unloadNamespace("robscale")

cat(sprintf("Platform: %s\n", Sys.info()["sysname"]))
cat(sprintf("Baseline L2/core: %d, has pdq_median_threshold: %s\n",
            cfg_baseline$l2_per_core,
            "pdq_median_threshold" %in% names(cfg_baseline)))
cat(sprintf("Adaptive L2/core: %d, pdq_median_threshold: %s\n",
            cfg_adaptive$l2_per_core,
            cfg_adaptive$pdq_median_threshold))
cat("\n")

K <- 1.482602218505602

# Sizes: regression range (3K-10K), crossover (~15K), speedup range (50K+)
SIZES <- c(1000, 3000, 5000, 10000, 15000, 20000, 50000, 100000, 500000)

results <- bench::press(
  n = SIZES,
  {
    set.seed(42 + n)
    x <- rnorm(n)
    iters <- if (n <= 10000) 500L else if (n <= 100000) 100L else 20L

    bench::mark(
      baseline = mad_baseline(x, constant = K),
      adaptive = mad_adaptive(x, constant = K),
      min_iterations = iters,
      check = TRUE
    )
  }
)

res <- as.data.frame(results)
res$expression <- as.character(res$expression)

keep <- c("expression", "n", "min", "median", "n_itr", "n_gc")
keep <- intersect(keep, names(res))
res <- res[, keep]

for (col in c("min", "median")) {
  if (col %in% names(res))
    res[[col]] <- as.numeric(res[[col]], units = "ms")
}

# Compute ratio: adaptive / baseline (> 1 means adaptive is faster)
res$ratio <- NA_real_
for (nn in unique(res$n)) {
  mask <- res$n == nn
  t_base <- res$median[mask & res$expression == "baseline"]
  t_adap <- res$median[mask & res$expression == "adaptive"]
  if (length(t_base) == 1 && length(t_adap) == 1) {
    res$ratio[mask & res$expression == "adaptive"] <- t_base / t_adap
    res$ratio[mask & res$expression == "baseline"] <- 1.0
  }
}

platform <- if (Sys.info()["sysname"] == "Darwin") "macos" else "linux"
outfile <- sprintf("adaptive_mad_ab_%s_%s.csv",
                   format(Sys.time(), "%Y%m%d"), platform)
write.csv(res, outfile, row.names = FALSE)
cat(sprintf("Results saved to: %s\n\n", outfile))

cat("=== A/B: adaptive vs baseline 0.2.0 (median ms) ===\n\n")
cat(sprintf("  %8s  %10s  %10s  %8s\n", "n", "baseline", "adaptive", "ratio"))
cat(paste(rep("-", 46), collapse = ""), "\n")
for (nn in SIZES) {
  t_b <- res$median[res$n == nn & res$expression == "baseline"]
  t_a <- res$median[res$n == nn & res$expression == "adaptive"]
  r   <- res$ratio[res$n == nn & res$expression == "adaptive"]
  cat(sprintf("  %8d  %10.4f  %10.4f  %8.3fx\n", nn, t_b, t_a, r))
}

# Gate: adaptive should not regress vs baseline at any size
adap_rows <- res[res$expression == "adaptive", ]
worst <- min(adap_rows$ratio, na.rm = TRUE)
worst_n <- adap_rows$n[which.min(adap_rows$ratio)]
cat(sprintf("\nWorst adaptive/baseline ratio: %.3fx at n=%d\n", worst, worst_n))
if (worst < 0.95) {
  cat("FAIL: regression detected (< 0.95x)\n")
} else if (worst < 1.00) {
  cat("WARN: marginal (< 1.00x but >= 0.95x)\n")
} else {
  cat("PASS: no regression vs baseline\n")
}
