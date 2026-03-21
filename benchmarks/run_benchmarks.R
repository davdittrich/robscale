# ---------------------------------------------------------------------------
# Build-time metadata helpers
# ---------------------------------------------------------------------------

#' Read COMMON_FLAGS / CXXFLAGS from the user's ~/.R/Makevars
#' Returns the raw flag string, or "" if not found.
read_user_cxxflags <- function() {
  mf <- file.path(path.expand("~"), ".R", "Makevars")
  if (!file.exists(mf)) return("")
  lines <- readLines(mf, warn = FALSE)
  # Strip comments and collapse line continuations
  lines <- sub("#.*", "", lines)
  lines <- lines[nchar(trimws(lines)) > 0]
  # Prefer COMMON_FLAGS (as used in this project's Makevars)
  hit <- grep("^COMMON_FLAGS\\s*=", lines, value = TRUE)
  if (!length(hit)) hit <- grep("^CXXFLAGS\\s*=", lines, value = TRUE)
  if (!length(hit)) return("")
  trimws(sub("^[A-Z_]+\\s*=\\s*", "", hit[1]))
}

#' Parse library detection from R CMD INSTALL log; combine with flags from
#' the configure-generated src/Makevars captured before source deletion.
#' @param install_log  Character vector: stdout+stderr from system2() install
#' @param pkg_cxxflags PKG_CXXFLAGS value read from tmp_src/src/Makevars
#' @return Named list: cxx_flags, tanh_backend, tbb_source, omp_simd (logical)
parse_build_info <- function(install_log, pkg_cxxflags = "") {
  log_str <- paste(install_log, collapse = "\n")

  # Full CXXFLAGS = user base flags + package-specific flags from configure
  user_flags <- read_user_cxxflags()
  parts <- trimws(c(user_flags, pkg_cxxflags))
  cxx_flags <- paste(parts[nchar(parts) > 0], collapse = " ")
  if (!nchar(cxx_flags)) cxx_flags <- NA_character_

  tanh_backend <- if (grepl("_ZGVdN4v_tanh detected|libmvec.*detected", log_str)) {
    "glibc libmvec (_ZGVdN4v_tanh)"
  } else if (grepl("SLEEF detected", log_str, fixed = TRUE)) {
    "SLEEF"
  } else if (grepl("Accelerate framework", log_str, fixed = TRUE)) {
    "Apple Accelerate (vvtanh)"
  } else if (grepl("fopenmp-simd supported", log_str, fixed = TRUE)) {
    "OpenMP SIMD"
  } else {
    "scalar std::tanh"
  }

  tbb_source <- if (grepl("System oneTBB detected via pkg-config", log_str, fixed = TRUE)) {
    "system oneTBB (pkg-config)"
  } else if (grepl("System oneTBB detected", log_str, fixed = TRUE)) {
    "system oneTBB (.so)"
  } else if (grepl("RcppParallel TBB detected", log_str, fixed = TRUE)) {
    "RcppParallel TBB"
  } else {
    "none (OpenMP parallel fallback)"
  }

  list(
    cxx_flags    = cxx_flags,
    tanh_backend = tanh_backend,
    tbb_source   = tbb_source,
    omp_simd     = grepl("fopenmp-simd supported", log_str, fixed = TRUE)
  )
}

#' Read CPU model name from /proc/cpuinfo or sysctl (macOS)
get_cpu_name <- function() {
  if (file.exists("/proc/cpuinfo")) {
    line <- system("grep 'model name' /proc/cpuinfo | head -n1 | cut -d':' -f2",
                   intern = TRUE)
    if (length(line)) return(trimws(line))
  }
  tryCatch(
    trimws(system("sysctl -n machdep.cpu.brand_string", intern = TRUE)),
    error = function(e) "Unknown CPU"
  )
}

#' Read scaling governor for cpu0 (Linux only)
get_cpu_governor <- function() {
  gov <- "/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"
  if (file.exists(gov)) trimws(readLines(gov, n = 1)) else "unknown"
}

# ---------------------------------------------------------------------------

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
  
  if (!is.null(attr(res_install, "status"))) {
    unlink(tmp_src, recursive = TRUE)
    stop("Installation failed: ", paste(res_install, collapse = "\n"))
  }

  # Read configure-generated src/Makevars BEFORE deleting the temp source.
  # This is the authoritative record of package-specific compile flags.
  pkg_makevars_path <- file.path(tmp_src, "src", "Makevars")
  pkg_cxxflags <- ""
  if (file.exists(pkg_makevars_path)) {
    mv_lines <- readLines(pkg_makevars_path, warn = FALSE)
    hit <- grep("^PKG_CXXFLAGS\\s*=", mv_lines, value = TRUE)
    if (length(hit)) pkg_cxxflags <- trimws(sub("^PKG_CXXFLAGS\\s*=\\s*", "", hit[1]))
  }

  unlink(tmp_src, recursive = TRUE)

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
  
  # Attach all build-time and system metadata captured at benchmark time
  res_obj$has_sleef <- has_sleef
  res_obj$build_info <- c(
    parse_build_info(res_install, pkg_cxxflags),
    list(
      pkg_version    = as.character(read.dcf(file.path(source_path, "DESCRIPTION"))[, "Version"]),
      r_version      = res_obj$sys_info$platform$version,
      platform       = res_obj$sys_info$platform$os,
      cpu_name       = get_cpu_name(),
      cpu_governor   = get_cpu_governor(),
      benchmark_date = as.character(Sys.Date())
    )
  )
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
    
    fast_mad <- function(x, constant = 1.4826) {
      m <- collapse::fmedian(x)
      return(collapse::fmedian(abs(x - m)) * constant)
    }

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
          mad_collapse = fast_mad(x),
          check = FALSE,
          min_iterations = get_min_iters(n),
          min_time = 1.0
        )
      }
    )

    list(
      revss = results_revss,
      robustbase = results_robustbase,
      new_estimators = results_new,
      pkg_versions = list(
        robustbase   = as.character(packageVersion("robustbase")),
        revss        = as.character(packageVersion("revss")),
        Hmisc        = as.character(packageVersion("Hmisc")),
        GiniDistance = as.character(packageVersion("GiniDistance")),
        collapse     = as.character(packageVersion("collapse"))
      )
    )
  }, args = list(lp = .libPaths()), show = TRUE)
  
  return(res_obj)
}
