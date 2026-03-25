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

#' Pre-flight system noise check (BM-4)
#' Measures the coefficient of variation (CV) of a minimal timing operation.
#' A CV > 10% indicates the system is too noisy for reliable micro-benchmarks.
#' @return Named list: cv (numeric), noisy (logical), warning (character or NULL)
check_system_noise <- function(reps = 1000L) {
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()[["elapsed"]]
    1L + 1L
    times[i] <- proc.time()[["elapsed"]] - t0
  }
  # Drop exact zeros (timer resolution) before computing CV
  times <- times[times > 0]
  cv <- if (length(times) > 1L) sd(times) / mean(times) else NA_real_
  noisy <- !is.na(cv) && cv > 0.10
  list(
    cv     = cv,
    noisy  = noisy,
    warning = if (noisy) sprintf(
      "System noise CV=%.1f%% > 10%% — results may be unreliable. Consider: performance governor, isolcpus, taskset.",
      cv * 100
    ) else NULL
  )
}

#' Pool time vectors from multiple bench::press runs (one per seed). (BM-2)
#' Returns the first result with the `time` list-column replaced by the
#' concatenated per-iteration times (in seconds, bench_time unit) from all
#' seeds, and all summary stats recomputed from the pooled distribution.
#' @param seed_results List of bench_mark/bench_press objects (one per seed).
pool_bench_press <- function(seed_results) {
  if (length(seed_results) == 1L) return(seed_results[[1L]])
  base <- seed_results[[1L]]
  # as.numeric() on bench_time returns seconds — pool across seeds per row.
  for (i in seq_len(nrow(base))) {
    pooled_s <- unlist(lapply(seed_results, function(sr) as.numeric(sr$time[[i]])))
    base$time[[i]] <- pooled_s   # plain numeric in seconds; as.numeric() in
  }                              # analyze_results.R handles it transparently.
  pool_s <- lapply(base$time, as.numeric)
  base$min      <- bench::as_bench_time(vapply(pool_s, min,    numeric(1L)))
  base$median   <- bench::as_bench_time(vapply(pool_s, median, numeric(1L)))
  base$mean     <- bench::as_bench_time(vapply(pool_s, mean,   numeric(1L)))
  base$max      <- bench::as_bench_time(vapply(pool_s, max,    numeric(1L)))
  # Update iteration counts and throughput from the pooled data.
  base$n_itr    <- vapply(pool_s, length, integer(1L))
  base[["itr/sec"]] <- 1 / vapply(pool_s, mean, numeric(1L))
  base$total_time <- bench::as_bench_time(vapply(pool_s, sum, numeric(1L)))
  base
}

# Seeds used for multi-seed benchmark runs (BM-2).
# Spaced 50 apart to avoid correlated .Random.seed states.
# 7 seeds (up from 3) to reduce variance in the pooled-median estimate at
# small n (n<=16), where the Aitken path hit vs miss ratio was biasing results.
BENCH_SEEDS <- c(42L, 92L, 142L, 192L, 242L, 292L, 342L)

# ---------------------------------------------------------------------------

#' Benchmark a specific version of robscale
#' @param install_env Named list of environment variables for installation
#' @param lib_path Path to a temporary library for installation
benchmark_robscale <- function(install_env = list(), lib_path = tempfile("lib_"), source_path = ".") {
  dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)

  # ── Pre-flight system checks (BM-4) ───────────────────────────────────────
  gov <- get_cpu_governor()
  if (!gov %in% c("performance", "unknown")) {
    message(sprintf(
      "WARNING: CPU governor is '%s', not 'performance'. Set with:\n  echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor",
      gov
    ))
  }
  noise <- check_system_noise()
  if (!is.null(noise$warning)) message("WARNING: ", noise$warning)

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
  res_obj <- callr::r(function(lib_path, lp, seeds) {
    # Isolate library paths for this process
    withr::with_libpaths(new = lib_path, action = "prefix", {
      library(robscale)
      library(bench)

      n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                   12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)
      n_large <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                   12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)

      get_min_iters <- function(n) {
        if (n <= 128) 10000L
        else if (n <= 2048) 2000L
        else if (n <= 16384) 500L
        else if (n <= 1048576) 20L
        else 5L
      }

      # Pool time vectors from multiple bench::press runs (one per seed). (BM-2)
      pool_bench_press <- function(seed_results) {
        if (length(seed_results) == 1L) return(seed_results[[1L]])
        base <- seed_results[[1L]]
        for (i in seq_len(nrow(base))) {
          pooled_s <- unlist(lapply(seed_results, function(sr) as.numeric(sr$time[[i]])))
          base$time[[i]] <- pooled_s
        }
        pool_s <- lapply(base$time, as.numeric)
        base$min         <- bench::as_bench_time(vapply(pool_s, min,    numeric(1L)))
        base$median      <- bench::as_bench_time(vapply(pool_s, median, numeric(1L)))
        base$mean        <- bench::as_bench_time(vapply(pool_s, mean,   numeric(1L)))
        base$max         <- bench::as_bench_time(vapply(pool_s, max,    numeric(1L)))
        base$n_itr       <- vapply(pool_s, length, integer(1L))
        base[["itr/sec"]]  <- 1 / vapply(pool_s, mean, numeric(1L))
        base$total_time  <- bench::as_bench_time(vapply(pool_s, sum, numeric(1L)))
        base
      }

      # ── Warmup (BM-3): stabilise CPU P-states and instruction cache ────────
      .wu_x <- rnorm(64L)
      for (.wu_i in seq_len(300L)) {
        robscale::robLoc(.wu_x)
        robscale::robScale(.wu_x)
        robscale::adm(.wu_x)
        robscale::qn(.wu_x)
        robscale::sn(.wu_x)
        robscale::gmd(.wu_x)
      }
      rm(.wu_x, .wu_i)
      gc(full = TRUE)

      # ── Small sample benchmarks (3-seed pooled) ────────────────────────────
      seed_results_m <- lapply(seeds, function(seed) {
        bench::press(
          n = n_small,
          {
            set.seed(seed + n)
            x <- rnorm(n)
            gc(full = TRUE)   # BM-3: flush GC before each configuration
            bench::mark(
              robLoc   = robscale::robLoc(x),
              robScale = robscale::robScale(x),
              adm      = robscale::adm(x),
              check    = FALSE,
              min_iterations = get_min_iters(n),
              min_time = 1.0
            )
          }
        )
      })
      results_m <- pool_bench_press(seed_results_m)

      # ── Large sample benchmarks (3-seed pooled) ────────────────────────────
      seed_results_scale <- lapply(seeds, function(seed) {
        bench::press(
          n = n_large,
          {
            set.seed(seed + n)
            x <- rnorm(n)
            gc(full = TRUE)
            bench::mark(
              qn = robscale::qn(x),
              sn = robscale::sn(x),
              check    = FALSE,
              min_iterations = get_min_iters(n),
              min_time = 1.0
            )
          }
        )
      })
      results_scale <- pool_bench_press(seed_results_scale)

      # ── New estimators (3-seed pooled) ─────────────────────────────────────
      seed_results_new <- lapply(seeds, function(seed) {
        bench::press(
          n = n_large,
          {
            set.seed(seed + n)
            x <- rnorm(n)
            gc(full = TRUE)
            bench::mark(
              gmd        = robscale::gmd(x),
              iqr_scaled = robscale::iqr_scaled(x),
              mad_scaled = robscale::mad_scaled(x),
              check    = FALSE,
              min_iterations = get_min_iters(n),
              min_time = 1.0
            )
          }
        )
      })
      results_new <- pool_bench_press(seed_results_new)

      # ── Ensemble estimator (scale_robust, auto_switch=FALSE) ───────────────
      # Benchmark the bootstrap ensemble path directly, bypassing the gmd()
      # fallback.  Ensemble is O(n_boot * n_estimators * estimator_cost), so
      # we cap at n=512 and use a lower iteration budget to keep runtime sane.
      # At n > 512, ensemble overhead would dominate the figure unproductively.
      n_ensemble <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512)
      get_min_iters_ens <- function(n) {
        if (n <=  16L) 1000L
        else if (n <= 64L)  200L
        else if (n <= 256L)  50L
        else 10L
      }
      seed_results_ensemble <- lapply(seeds, function(seed) {
        bench::press(
          n = n_ensemble,
          {
            set.seed(seed + n)
            x <- rnorm(n)
            gc(full = TRUE)
            bench::mark(
              scale_robust = robscale::scale_robust(x, auto_switch = FALSE),
              check    = FALSE,
              min_iterations = get_min_iters_ens(n),
              min_time = 1.0
            )
          }
        )
      })
      results_ensemble <- pool_bench_press(seed_results_ensemble)

      list(
        m_estimators     = results_m,
        scale_estimators = results_scale,
        new_estimators   = results_new,
        ensemble         = results_ensemble,
        sys_info         = sessioninfo::session_info()
      )
    })
  }, args = list(lib_path = lib_path, lp = .libPaths(), seeds = BENCH_SEEDS), show = TRUE)
  
  # Attach all build-time and system metadata captured at benchmark time
  res_obj$has_sleef <- has_sleef
  res_obj$build_info <- c(
    parse_build_info(res_install, pkg_cxxflags),
    list(
      pkg_version    = as.character(read.dcf(file.path(source_path, "DESCRIPTION"))[, "Version"]),
      r_version      = res_obj$sys_info$platform$version,
      platform       = res_obj$sys_info$platform$os,
      cpu_name       = get_cpu_name(),
      cpu_governor   = gov,             # captured in pre-flight check above
      noise_cv       = noise$cv,        # system noise level at benchmark time
      bench_seeds    = BENCH_SEEDS,     # seeds used for multi-seed pooling
      benchmark_date = as.character(Sys.Date())
    )
  )
  return(res_obj)
}

#' Benchmark legacy packages
benchmark_legacy <- function() {
  # Also run legacy in its own process for consistency
  res_obj <- callr::r(function(lp, seeds) {
    library(bench)
    library(robustbase)
    library(revss)
    library(Hmisc)
    library(GiniDistance)
    library(collapse)

    n_small <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                 12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)
    n_large <- c(3, 4, 5, 6, 7, 8, 10, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
                 12288, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 10000000)

    get_min_iters <- function(n) {
      if (n <= 128) 10000L
      else if (n <= 2048) 2000L
      else if (n <= 16384) 500L
      else if (n <= 1048576) 20L
      else 5L
    }

    # Pool time vectors from multiple bench::press runs (one per seed). (BM-2)
    pool_bench_press <- function(seed_results) {
      if (length(seed_results) == 1L) return(seed_results[[1L]])
      base <- seed_results[[1L]]
      for (i in seq_len(nrow(base))) {
        pooled_ns <- unlist(lapply(seed_results, function(sr) as.numeric(sr$time[[i]])))
        base$time[[i]] <- pooled_ns
      }
      pool_ns <- lapply(base$time, as.numeric)
      base$min    <- bench::as_bench_time(vapply(pool_ns, min,    numeric(1L)))
      base$median <- bench::as_bench_time(vapply(pool_ns, median, numeric(1L)))
      base$mean   <- bench::as_bench_time(vapply(pool_ns, mean,   numeric(1L)))
      base$max    <- bench::as_bench_time(vapply(pool_ns, max,    numeric(1L)))
      base
    }

    # ── Warmup (BM-3) ──────────────────────────────────────────────────────
    .wu_x <- rnorm(64L)
    for (.wu_i in seq_len(300L)) {
      revss::robLoc(.wu_x)
      revss::robScale(.wu_x)
      robustbase::Qn(.wu_x)
      robustbase::Sn(.wu_x)
    }
    rm(.wu_x, .wu_i)
    gc(full = TRUE)

    fast_mad <- function(x, constant = 1.4826) {
      m <- collapse::fmedian(x)
      collapse::fmedian(abs(x - m)) * constant
    }

    # ── revss (3-seed pooled) ──────────────────────────────────────────────
    seed_results_revss <- lapply(seeds, function(seed) {
      bench::press(
        n = n_small,
        {
          set.seed(seed + n)
          x <- rnorm(n)
          gc(full = TRUE)
          bench::mark(
            robLoc   = revss::robLoc(x),
            robScale = revss::robScale(x),
            adm      = revss::adm(x),
            check    = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0
          )
        }
      )
    })
    results_revss <- pool_bench_press(seed_results_revss)

    # ── robustbase (3-seed pooled) ─────────────────────────────────────────
    seed_results_robustbase <- lapply(seeds, function(seed) {
      bench::press(
        n = n_large,
        {
          set.seed(seed + n)
          x <- rnorm(n)
          gc(full = TRUE)
          bench::mark(
            qn = robustbase::Qn(x),
            sn = robustbase::Sn(x),
            check    = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0
          )
        }
      )
    })
    results_robustbase <- pool_bench_press(seed_results_robustbase)

    # ── New estimator competitors (3-seed pooled) ──────────────────────────
    seed_results_new <- lapply(seeds, function(seed) {
      bench::press(
        n = n_large,
        {
          set.seed(seed + n)
          x <- rnorm(n)
          gc(full = TRUE)
          bench::mark(
            gmd          = Hmisc::GiniMd(x) * 0.886226925452758,
            gmd_gd       = GiniDistance::gmd(x) * 0.886226925452758,
            iqr_scaled   = stats::IQR(x) * 0.741301109252801,
            iqr_collapse = diff(collapse::fquantile(x, c(0.25, 0.75))) * 0.741301109252801,
            mad_scaled   = stats::mad(x),
            mad_collapse = fast_mad(x),
            check    = FALSE,
            min_iterations = get_min_iters(n),
            min_time = 1.0
          )
        }
      )
    })
    results_new <- pool_bench_press(seed_results_new)

    list(
      revss          = results_revss,
      robustbase     = results_robustbase,
      new_estimators = results_new,
      pkg_versions   = list(
        robustbase   = as.character(packageVersion("robustbase")),
        revss        = as.character(packageVersion("revss")),
        Hmisc        = as.character(packageVersion("Hmisc")),
        GiniDistance = as.character(packageVersion("GiniDistance")),
        collapse     = as.character(packageVersion("collapse"))
      )
    )
  }, args = list(lp = .libPaths(), seeds = BENCH_SEEDS), show = TRUE)

  return(res_obj)
}
