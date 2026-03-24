## Pre-RS0 baseline capture for robScale()
## Run BEFORE any WU-RS0 code changes. Uses robScale() directly.
## Output: bench/rob_scale_pre_rs0_baseline.rds

library(robscale)
library(bench)

sizes_pre <- c(4, 5, 8, 10, 15, 16, 17, 30, 50, 64, 65, 100, 200, 500, 1000,
               4096, 8192, 32768, 65536)

pre_results <- lapply(sizes_pre, function(n) {
  set.seed(42)
  x <- rnorm(n)
  bm <- bench::mark(
    robScale(x),
    min_iterations = if (n <= 100) 500L else if (n <= 1000) 200L else 50L,
    max_time = 30,
    check = FALSE
  )
  list(n = n, median_ns = as.numeric(bm$median) * 1e9)
})

saveRDS(pre_results, file.path("bench", "rob_scale_pre_rs0_baseline.rds"))

cat("Pre-RS0 baseline captured:\n")
for (r in pre_results) {
  cat(sprintf("  n = %6d  median = %.1f ns\n", r$n, r$median_ns))
}
cat("Saved to bench/rob_scale_pre_rs0_baseline.rds\n")
