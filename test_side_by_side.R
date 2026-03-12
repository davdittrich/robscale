library(bench)

# Default is optimized in main library
fast_qn <- robscale::qn

# Create a slow library
lib_slow <- tempdir()
dir.create(lib_slow, showWarnings = FALSE)

Sys.setenv(ROBSCALE_FAST = "0")
install.packages(".", lib = lib_slow, repos = NULL, type = "source", INSTALL_opts = "--no-multiarch", quiet = TRUE)
Sys.unsetenv("ROBSCALE_FAST")

# Load slow package into a new namespace temporarily to avoid conflict
slow_qn <- callr::r(function(lib, n) {
  .libPaths(c(lib, .libPaths()))
  x <- rnorm(n)
  res <- bench::mark(
    qn_slow = robscale::qn(x),
    min_iterations = 1000, 
    check = FALSE
  )
  as.numeric(res$median)
}, args = list(lib = lib_slow, n = 3))

fast_time <- as.numeric(bench::mark(qn_fast = fast_qn(rnorm(3)), min_iterations = 1000, check = FALSE)$median)

cat("Slow median:", slow_qn * 1e6, "us\n")
cat("Fast median:", fast_time * 1e6, "us\n")
