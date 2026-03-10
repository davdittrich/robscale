#' Benchmark a specific version of robscale
#' @param install_env Named list of environment variables for installation
#' @param lib_path Path to a temporary library for installation
benchmark_robscale <- function(install_env = list(), lib_path = tempfile("lib_")) {
  dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)
  
  # Install the package into the temporary library
  # We need to ensure the installer sees the current library paths (renv)
  current_libs <- .libPaths()
  message("Current libPaths: ", paste(current_libs, collapse = ", "))
  
  # Set R_LIBS for the R CMD INSTALL process
  env_libs <- paste(current_libs, collapse = .Platform$path.sep)
  envs <- c(install_env, R_LIBS = env_libs, R_LIBS_USER = env_libs)
  
  # Using R CMD INSTALL directly for better control
  res <- system2(R.home("bin/R"), 
                args = c("CMD", "INSTALL", ".", 
                        paste0("--library=", lib_path),
                        "--no-multiarch", "--no-test-load", "--quiet"),
                env = paste0(names(envs), "=", envs),
                stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(res, "status"))) {
    stop("Installation failed: ", paste(res, collapse = "\n"))
  }
  
  # Load the library and run benchmarks
  # IMPORTANT: action = "prefix" ensures we keep the original (renv) libraries where dependencies are!
  withr::with_libpaths(new = lib_path, action = "prefix", {
    library(robscale)
    library(bench)
    
    # Define sizes
    n_small <- c(5, 10, 20)
    n_large <- c(100, 1000, 10000, 100000)
    
    # Small sample benchmarks (vs revss if installed, or just internal)
    results_m <- bench::press(
      n = n_small,
      {
        x <- rnorm(n)
        bench::mark(
          robLoc = robscale::robLoc(x),
          robScale = robscale::robScale(x),
          adm = robscale::adm(x),
          check = FALSE,
          min_iterations = 100
        )
      }
    )
    
    # Large sample benchmarks
    results_scale <- bench::press(
      n = n_large,
      {
        x <- rnorm(n)
        bench::mark(
          qn = robscale::qn(x),
          sn = robscale::sn(x),
          check = FALSE,
          min_iterations = 20
        )
      }
    )
    
    list(
      m_estimators = results_m,
      scale_estimators = results_scale,
      sys_info = sessioninfo::session_info()
    )
  })
}

#' Benchmark legacy packages
benchmark_legacy <- function() {
  library(bench)
  library(robustbase)
  library(revss)
  
  n_small <- c(5, 10, 20)
  n_large <- c(100, 1000, 10000, 100000)
  
  results_revss <- bench::press(
    n = n_small,
    {
      x <- rnorm(n)
      bench::mark(
        robLoc = revss::robLoc(x),
        robScale = revss::robScale(x),
        adm = revss::adm(x),
        check = FALSE,
        min_iterations = 100
      )
    }
  )
  
  results_robustbase <- bench::press(
    n = n_large,
    {
      x <- rnorm(n)
      bench::mark(
        qn = robustbase::Qn(x),
        sn = robustbase::Sn(x),
        check = FALSE,
        min_iterations = 20
      )
    }
  )
  
  list(
    revss = results_revss,
    robustbase = results_robustbase
  )
}
