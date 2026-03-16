#!/usr/bin/env Rscript --vanilla
library(robscale)

cat("=== Sort Shootout: n <= 64 ===\n")
cat("Platform:", Sys.info()["sysname"], Sys.info()["machine"], "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -- Correctness --
cat("Checking correctness... ")
corr <- robscale:::sort_shootout_correctness()
failures <- corr[!corr$matches_reference, ]
if (nrow(failures) > 0) {
  cat("FAILED\n")
  print(failures)
  stop("Correctness check failed!")
}
cat("OK (", nrow(corr), " checks passed)\n\n", sep = "")

# -- Benchmark --
cat("Running benchmark (1000 x 101 rounds per cell)...\n")
bm <- robscale:::sort_shootout_benchmark(batch_size = 1000L, nrounds = 101L)

# Pivot to wide format: rows = n, columns = method
wide <- reshape(bm[, c("method", "n", "median_ns")],
                direction = "wide", idvar = "n", timevar = "method",
                v.names = "median_ns")
names(wide) <- sub("^median_ns\\.", "", names(wide))

# Find winner per n
methods <- setdiff(names(wide), "n")
wide$winner <- methods[apply(wide[methods], 1, which.min)]

# Speedup of winner vs std_sort
wide$speedup <- sapply(seq_len(nrow(wide)), function(i) {
  round(wide$std_sort[i] / wide[[wide$winner[i]]][i], 2)
})

cat("\nMedian time per sort (ns):\n")
fmt <- wide
for (m in methods) fmt[[m]] <- sprintf("%7.1f", fmt[[m]])
print(fmt, row.names = FALSE)

# Save
outfile <- paste0("sort_shootout_", format(Sys.Date(), "%Y%m%d"), ".csv")
write.csv(bm, outfile, row.names = FALSE)
cat("\nFull results saved to", outfile, "\n")
