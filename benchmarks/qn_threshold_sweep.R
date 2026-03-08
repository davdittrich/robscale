library(bench)
library(robscale)
library(fastqnsn)
library(ggplot2)

# Grid of samples sizes around the crossover
n_grid <- c(256, 384, 448, 512, 640, 768, 896, 1024, 1536, 2048)

# Get current threshold from the package
config <- robscale:::get_qnsn_config()
current_threshold <- config$qn_exact_threshold
cat(sprintf("Benchmarking with current QN_EXACT_THRESHOLD: %d\n", current_threshold))

results <- bench::press(n = n_grid, {
  set.seed(42)
  x <- rnorm(n)
  bench::mark(
    robscale = robscale::qn(x),
    fastqnsn = fastqnsn::qn(x),
    min_iterations = 200,
    check = TRUE
  )
})

# Save results to a CSV named after the threshold
res_flat <- as.data.frame(results)
res_flat$expression <- as.character(results$expression)
res_flat$median <- as.numeric(results$median)
res_flat$min <- as.numeric(results$min)
res_flat$itr_sec <- as.numeric(results[["itr/sec"]])

# Remove list columns that cause write.csv to fail
res_flat <- res_flat[, !sapply(res_flat, is.list)]

out_file <- sprintf("benchmarks/qn_threshold_sweep_%d.csv", current_threshold)
write.csv(res_flat, out_file, row.names = FALSE)

cat(sprintf("Results saved to %s\n", out_file))

# Simple comparison plot
p <- ggplot(results, aes(x = n, y = as.numeric(median), color = expression)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(title = sprintf("Qn Performance (Threshold = %d)", current_threshold),
       y = "Median Time (sec)", x = "Sample Size (n)") +
  theme_minimal()

ggsave(sprintf("benchmarks/qn_threshold_sweep_%d.png", current_threshold), p, width = 8, height = 6)
