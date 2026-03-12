library(bench)
fast_adm <- robscale::adm

# Load slow package into a new namespace temporarily to avoid conflict
lib_slow <- tempdir()
Sys.setenv(ROBSCALE_FAST = "0")
install.packages(".", lib = lib_slow, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch", quiet = TRUE)
Sys.unsetenv("ROBSCALE_FAST")
slow_adm <- callr::r(function(lib) {
  .libPaths(c(lib, .libPaths()))
  function(x) robscale::adm(x)
}, args = list(lib = lib_slow))

set.seed(42)
x <- rnorm(256)

fast_time <- as.numeric(bench::mark(fast_adm(x), min_iterations=1000, check=FALSE)$median)
slow_time <- as.numeric(bench::mark(slow_adm(x), min_iterations=1000, check=FALSE)$median)

cat("Slow median:", slow_time * 1e6, "us\n")
cat("Fast median:", fast_time * 1e6, "us\n")
