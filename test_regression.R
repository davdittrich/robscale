library(bench)
library(robscale)

res <- bench::press(
  n = c(3, 4, 10, 64, 128),
  {
    set.seed(n)
    x <- rnorm(n)
    bench::mark(
      qn_robscale = robscale::qn(x),
      min_iterations = 1000,
      check = FALSE
    )
  }
)
print(res)
