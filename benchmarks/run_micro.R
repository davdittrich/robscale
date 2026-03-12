source("benchmarks/run_benchmarks.R")

run_micro <- function() {
  cat("Building fast_mode for micro benchmark...\n")
  fast_lib <- tempfile("lib_fast_micro_")
  dir.create(fast_lib, recursive = TRUE, showWarnings = FALSE)
  
  # Install
  envs <- c(list(ROBSCALE_FAST=1), R_LIBS = .libPaths()[1], R_LIBS_USER = .libPaths()[1])
  res_install <- system2(R.home("bin/R"), 
                args = c("CMD", "INSTALL", ".", 
                        paste0("--library=", fast_lib),
                        "--no-multiarch", "--no-test-load"),
                env = paste0(names(envs), "=", envs),
                stdout = TRUE, stderr = TRUE)
                
  callr::r(function(lib_path, lp) {
    withr::with_libpaths(new = lib_path, action = "prefix", {
      library(robscale)
      library(revss)
      library(bench)
      
      cat("\n--- adm n=4 ---\n")
      set.seed(42 + 4)
      x <- rnorm(4)
      res <- bench::mark(
        revss = revss::adm(x),
        robscale_fast = robscale::adm(x),
        check = FALSE, min_iterations = 10000, min_time = 0.5
      )
      print(res)
      cat(sprintf("Speedup: %.3f\n", as.numeric(res$median[1]) / as.numeric(res$median[2])))
      
      cat("\n--- robLoc n=5 ---\n")
      set.seed(42 + 5)
      x <- rnorm(5)
      res <- bench::mark(
        revss = revss::robLoc(x),
        robscale_fast = robscale::robLoc(x),
        check = FALSE, min_iterations = 10000, min_time = 0.5
      )
      print(res)
      cat(sprintf("Speedup: %.3f\n", as.numeric(res$median[1]) / as.numeric(res$median[2])))

      cat("\n--- robLoc n=6 ---\n")
      set.seed(42 + 6)
      x <- rnorm(6)
      res <- bench::mark(
        revss = revss::robLoc(x),
        robscale_fast = robscale::robLoc(x),
        check = FALSE, min_iterations = 10000, min_time = 0.5
      )
      print(res)
      cat(sprintf("Speedup: %.3f\n", as.numeric(res$median[1]) / as.numeric(res$median[2])))
      
      cat("\n--- robScale n=5 ---\n")
      set.seed(42 + 5)
      x <- rnorm(5)
      res <- bench::mark(
        revss = revss::robScale(x),
        robscale_fast = robscale::robScale(x),
        check = FALSE, min_iterations = 10000, min_time = 0.5
      )
      print(res)
      cat(sprintf("Speedup: %.3f\n", as.numeric(res$median[1]) / as.numeric(res$median[2])))
    })
  }, args = list(lib_path = fast_lib, lp = .libPaths()), show = TRUE)
}

run_micro()
