library(bench)
fast_qn <- robscale::qn

# Load slow package into a new namespace temporarily to avoid conflict
lib_slow <- tempdir()
Sys.setenv(ROBSCALE_FAST = "0")
install.packages(".", lib = lib_slow, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch", quiet = TRUE)
Sys.unsetenv("ROBSCALE_FAST")
slow_qn <- callr::r(function(lib) {
  .libPaths(c(lib, .libPaths()))
  function(x) robscale::qn(x)
}, args = list(lib = lib_slow))

set.seed(42)
x <- rnorm(100000)

res_fast <- bench::mark(fast_qn(x), min_iterations=10, check=FALSE)
res_slow <- bench::mark(slow_qn(x), min_iterations=10, check=FALSE)

cat("Slow median:", as.numeric(res_slow$median) * 1e3, "ms\n")
cat("Fast median:", as.numeric(res_fast$median) * 1e3, "ms\n")
