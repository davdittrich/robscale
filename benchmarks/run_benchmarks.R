#' Benchmark a specific version of robscale
#' @param install_env Named list of environment variables for installation
#' @param lib_path Path to a temporary library for installation
benchmark_robscale <- function(install_env = list(), lib_path = tempfile("lib_"), source_path = ".") {
  dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)

  # Create a temporary source directory to avoid polluting the current src/ 
  # or using stale .o files from a different build Type
  tmp_src <- tempfile("src_")
  dir.create(tmp_src, showWarnings = FALSE, recursive = TRUE)
  # Copy source_path to tmp_src
  file.copy(file.path(source_path, "DESCRIPTION"), tmp_src)
  file.copy(file.path(source_path, "NAMESPACE"), tmp_src, overwrite = TRUE)
  file.copy(file.path(source_path, "R"), tmp_src, recursive = TRUE)
  file.copy(file.path(source_path, "src"), tmp_src, recursive = TRUE)
  file.copy(file.path(source_path, "inst"), tmp_src, recursive = TRUE)
  file.copy(file.path(source_path, "configure"), tmp_src)
  file.copy(file.path(source_path, "cleanup"), tmp_src)
  
  # CRITICAL: Purge any object files or shared libraries that might have been 
  # copied from the workspace to the temporary source.
  unlink(list.files(file.path(tmp_src, "src"), 
                    pattern = "\\.o$|\\.so$|\\.dll$", 
                    full.names = TRUE))

  # Install the package into the temporary library from the temporary source
  current_libs <- .libPaths()
  env_libs <- paste(current_libs, collapse = .Platform$path.sep)
  envs <- c(install_env, R_LIBS = env_libs, R_LIBS_USER = env_libs)
  
  res_install <- system2(R.home("bin/R"), 
                args = c("CMD", "INSTALL", tmp_src, 
                        paste0("--library=", lib_path),
                        "--no-multiarch", "--no-test-load"),
                env = paste0(names(envs), "=", envs),
                stdout = TRUE, stderr = TRUE)
  
  unlink(tmp_src, recursive = TRUE)

  if (!is.null(attr(res_install, "status"))) {
    stop("Installation failed: ", paste(res_install, collapse = "\n"))
  }
  
  # Determine if SLEEF was used from the installation log
  has_sleef <- any(grepl("SLEEF detected|-DROBSCALE_HAS_SLEEF|-lsleef", res_install, ignore.case = TRUE))

  # Run benchmarks in a completely isolated R subprocess to prevent DLL/namespace 
  # caching across targets in the pipeline.
  res_obj <- callr::r(function(lib_path, lp) {
    # Isolate library paths for this process
    withr::with_libpaths(new = lib_path, action = "prefix", {
      library(robscale)
      library(bench)
      
      n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)
      n_large <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 
                   12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)
      
      get_min_iters <- function(n) {
        if (n <= 128) 10000L
        else if (n <= 2048) 2000L
        else if (n <= 16384) 500L
        else if (n <= 1048576) 20L
        else 5L
      }
      
      # Small sample benchmarks
      results_m <- bench::press(
        n = n_small,
        {
          set.seed(42 + n)
          x <- rnorm(n)
          bench::mark(
            robLoc = robscale::robLoc(x),
            robScale = robscale::robScale(x),
            adm = robscale::adm(x),
            check = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0 # Beat noise
          )
        }
      )
      
      # Large sample benchmarks
      results_scale <- bench::press(
        n = n_large,
        {
          set.seed(42 + n)
          x <- rnorm(n)
          bench::mark(
            qn = robscale::qn(x),
            sn = robscale::sn(x),
            check = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0
          )
        }
      )
      
      # New estimators benchmarks (gmd, iqr_scaled, mad_scaled)
      results_new <- bench::press(
        n = n_large,
        {
          set.seed(42 + n)
          x <- rnorm(n)
          bench::mark(
            gmd = robscale::gmd(x),
            iqr_scaled = robscale::iqr_scaled(x),
            mad_scaled = robscale::mad_scaled(x),
            check = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0
          )
        }
      )

      list(
        m_estimators = results_m,
        scale_estimators = results_scale,
        new_estimators = results_new,
        sys_info = sessioninfo::session_info()
      )
    })
  }, args = list(lib_path = lib_path, lp = .libPaths()), show = TRUE)
  
  res_obj$has_sleef <- has_sleef
  return(res_obj)
}

#' Benchmark legacy packages
benchmark_legacy <- function() {
  # Also run legacy in its own process for consistency
  res_obj <- callr::r(function(lp) {
    library(bench)
    library(robustbase)
    library(revss)
    library(Hmisc)
    library(GiniDistance)
    library(collapse)
    
    n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384)
    n_large <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 
                 12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)
    
    get_min_iters <- function(n) {
        if (n <= 128) 10000L
        else if (n <= 2048) 2000L
        else if (n <= 16384) 500L
        else if (n <= 1048576) 20L
        else 5L
    }
    
    results_revss <- bench::press(
      n = n_small,
      {
        set.seed(42 + n)
        x <- rnorm(n)
        bench::mark(
          robLoc = revss::robLoc(x),
          robScale = revss::robScale(x),
          adm = revss::adm(x),
          check = FALSE,
          min_iterations = get_min_iters(n),
          min_time = 1.0
        )
      }
    )
    
    results_robustbase <- bench::press(
      n = n_large,
      {
        set.seed(42 + n)
        x <- rnorm(n)
        bench::mark(
          qn = robustbase::Qn(x),
          sn = robustbase::Sn(x),
          check = FALSE,
          min_iterations = get_min_iters(n),
          min_time = 1.0
        )
      }
    )
    
    # New estimator competitors
    results_new <- bench::press(
      n = n_large,
      {
        set.seed(42 + n)
        x <- rnorm(n)
        bench::mark(
          gmd = Hmisc::GiniMd(x) * 0.886226925452758,
          gmd_gd = GiniDistance::gmd(x) * 0.886226925452758,
          iqr_scaled = stats::IQR(x) * 0.741301109252801,
          iqr_collapse = diff(collapse::fquantile(x, c(0.25, 0.75))) * 0.741301109252801,
          mad_scaled = stats::mad(x),
          mad_collapse = collapse::fmad(x, na.rm = FALSE),
          check = FALSE,
          min_iterations = get_min_iters(n),
          min_time = 1.0
        )
      }
    )

    list(
      revss = results_revss,
      robustbase = results_robustbase,
      new_estimators = results_new
    )
  }, args = list(lp = .libPaths()), show = TRUE)
  
  return(res_obj)
}
