# Threshold calibration: robScale median/MAD selection crossover (FR vs pdqselect)
#
# Purpose: measure the n at which pdqselect overtakes Floyd-Rivest for the
# median and MAD selection calls inside rob_scale.cpp.
#
# Working-set context: rob_scale uses 1 active array + ~1 warm array (x),
# giving lighter cache pressure than IQR/MAD (divisor ~2, threshold ~32-65K on Ryzen).
#
# Usage:
#   Rscript benchmarks/thresh_robscale_selection.R
#   # On remote: scp to 192.168.1.43 then ssh dd@192.168.1.43 'Rscript ~/benchmarks/thresh_robscale_selection.R'

suppressPackageStartupMessages({
  library(bench)
  library(robscale)
})

pkg_dir <- if (file.exists("DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(pkg_dir, "DESCRIPTION")))

# Install from current working tree before benchmarking
message("Installing current package...")
system2("R", c("CMD", "INSTALL", "--no-multiarch", "--no-test-load", pkg_dir),
        stdout = FALSE, stderr = FALSE)

sizes <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L, 100000L, 500000L)
n_iter <- 30L

distributions <- list(
  gaussian    = function(n) { set.seed(42); rnorm(n) },
  uniform     = function(n) { set.seed(42); runif(n) },
  contaminated = function(n) {
    set.seed(42)
    x <- rnorm(n)
    x[sample(n, max(1L, n %/% 10L))] <- rnorm(max(1L, n %/% 10L), mean = 10, sd = 3)
    x
  }
)

results <- vector("list", length(sizes) * length(distributions))
k <- 1L

for (dist_name in names(distributions)) {
  for (n in sizes) {
    x <- distributions[[dist_name]](n)
    bm <- bench::mark(
      robScale(x),
      min_iterations = n_iter,
      memory = FALSE
    )
    results[[k]] <- data.frame(
      distribution = dist_name,
      n = n,
      median_us = as.numeric(bm$median) * 1e6
    )
    k <- k + 1L
    cat(sprintf("  robScale  dist=%-12s  n=%7d  %.1f us\n",
                dist_name, n, results[[k - 1L]]$median_us))
  }
}

df <- do.call(rbind, results)

out_file <- sprintf("benchmarks/thresh_robscale_%s_%s.csv",
                    format(Sys.Date(), "%Y%m%d"),
                    ifelse(Sys.info()["sysname"] == "Darwin", "macos", "linux"))
write.csv(df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

# Report current threshold
cfg <- get_qnsn_config()
cat(sprintf("\nCurrent pdq_robscale_threshold: %d\n", cfg$pdq_robscale_threshold))
cat(sprintf("L2 per core: %d bytes\n", cfg$l2_per_core))
cat("\nNote: crossover n cannot be read directly from robScale timing alone\n")
cat("(both FR and pdqselect paths are invoked depending on n vs threshold).\n")
cat("To isolate paths, build two versions forcing FR-only vs pdqselect-only.\n")
