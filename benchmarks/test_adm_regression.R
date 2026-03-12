# benchmarks/test_adm_regression.R
library(bench)
library(ggplot2)

# Build and load the package with fast optimizations enabled
Sys.setenv(ROBSCALE_FAST = 1)
devtools::load_all(".", quiet = TRUE)

# We are testing the n=4 regression for adm
n <- 4
set.seed(42 + n)
x <- rnorm(n)

cat(sprintf("=== Benchmarking adm for n=%d ===\n", n))
res <- bench::mark(
  revss = revss::adm(x),
  robscale_fast = robscale::adm(x),
  check = FALSE,
  min_time = 2,
  min_iterations = 10000
)

print(res)
speedup <- as.numeric(res$median[1]) / as.numeric(res$median[2])
cat(sprintf("\nSpeedup (robscale_fast / revss): %.3fx\n", speedup))
if (speedup < 1.0) {
  cat("❌ REGRESSION DETECTED.\n")
} else {
  cat("✅ NO REGRESSION.\n")
}

# We also test robLoc at n=5 and n=6
for (n in c(5, 6)) {
  set.seed(42 + n)
  x <- rnorm(n)
  
  cat(sprintf("\n=== Benchmarking robLoc for n=%d ===\n", n))
  res <- bench::mark(
    revss = revss::robLoc(x),
    robscale_fast = robscale::robLoc(x),
    check = FALSE,
    min_time = 2,
    min_iterations = 10000
  )
  
  print(res)
  speedup <- as.numeric(res$median[1]) / as.numeric(res$median[2])
  cat(sprintf("\nSpeedup (robscale_fast / revss): %.3fx\n", speedup))
  if (speedup < 1.0) {
    cat("❌ REGRESSION DETECTED.\n")
  } else {
    cat("✅ NO REGRESSION.\n")
  }
}
