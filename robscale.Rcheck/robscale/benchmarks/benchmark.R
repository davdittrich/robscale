if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Please install microbenchmark: install.packages('microbenchmark')")
}
if (!requireNamespace("revss", quietly = TRUE)) {
  stop("Please install revss for comparison benchmarks")
}

library(microbenchmark)
library(revss)
library(robscale)

cat("=== robscale vs revss benchmark ===\n\n")

for (n in c(3, 4, 5, 6, 7, 8, 20, 100)) {
  set.seed(42)
  x <- runif(n, -100, 100)

  cat(sprintf("--- n = %d ---\n", n))

  mb_adm <- microbenchmark(
    revss = revss::adm(x),
    robscale = robscale::adm(x),
    times = 10000L
  )
  cat("adm:\n")
  print(summary(mb_adm)[, c("expr", "median", "mean")])

  mb_loc <- microbenchmark(
    revss = revss::robLoc(x),
    robscale = robscale::robLoc(x),
    times = 10000L
  )
  cat("robLoc:\n")
  print(summary(mb_loc)[, c("expr", "median", "mean")])

  mb_scale <- microbenchmark(
    revss = revss::robScale(x),
    robscale = robscale::robScale(x),
    times = 10000L
  )
  cat("robScale:\n")
  print(summary(mb_scale)[, c("expr", "median", "mean")])
  cat("\n")
}
