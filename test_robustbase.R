library(bench)
library(robustbase)
library(robscale)

x <- rnorm(3)
res <- bench::mark(
  qn_robscale = robscale::qn(x),
  qn_robustbase = robustbase::Qn(x),
  min_iterations = 5000,
  check = FALSE
)
print(res)
