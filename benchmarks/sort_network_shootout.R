#!/usr/bin/env Rscript
# Sort network shootout: benchmark Knuth vs alternatives for n<=8
#
# Usage: Rscript benchmarks/sort_network_shootout.R
# Requires: robscale installed

library(robscale)

cat("Running sort shootout benchmark (10000 reps per method x n) ...\n")
bm <- robscale:::sort_shootout_benchmark()

# Pivot to wide format: one column per method
wide <- reshape(bm[, c("method", "n", "median_ns")],
                idvar = "n", timevar = "method",
                direction = "wide")
names(wide) <- sub("median_ns\\.", "", names(wide))

cat("\n", strrep("=", 70), "\n")
cat("  Median time (ns) by method and n\n")
cat(strrep("=", 70), "\n")
print(wide, row.names = FALSE)

# Find winner per n
method_names <- unique(bm$method)
cat("\n", strrep("=", 70), "\n")
cat("  Winner per n\n")
cat(strrep("=", 70), "\n")
for (nn in sort(unique(bm$n))) {
  sub <- bm[bm$n == nn, ]
  best <- sub[which.min(sub$median_ns), ]
  knuth_time <- sub$median_ns[sub$method == "knuth_network"]
  speedup <- knuth_time / best$median_ns
  cat(sprintf("  n=%d: %-16s (%.1f ns, %.2fx vs knuth)\n",
              nn, best$method, best$median_ns, speedup))
}

# Save CSV
out_file <- file.path("benchmarks",
                      paste0("sort_shootout_", format(Sys.Date(), "%Y%m%d"), ".csv"))
write.csv(bm, out_file, row.names = FALSE)
cat(sprintf("\nResults saved to %s\n", out_file))
