#!/usr/bin/env Rscript
# Three-way benchmark: revss v2.0.0 (pure R) vs revss v3.0.0 (Fortran) vs
# robscale (C++/SIMD). Measures speed only — numerical differences are ignored.
#
# Usage: Rscript benchmarks/bench_revss_v2_v3.R
# Expected runtime: ~3-5 min

library(bench)

# ---------------------------------------------------------------------------
# Grid & iteration schedule (matches quick_regression_check.R)
# ---------------------------------------------------------------------------
sample_sizes <- c(4L, 8L, 16L, 64L, 256L, 1024L, 4096L, 16384L)

get_iters <- function(n) {
  if (n <= 128L)        200L
  else if (n <= 4096L)  100L
  else if (n <= 65536L)  50L
  else                    20L
}

# ---------------------------------------------------------------------------
# Install revss v2.0.0 from CRAN archive into a temp library
# ---------------------------------------------------------------------------
install_revss_v2 <- function(lib_path) {
  dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)
  url <- "https://cran.r-project.org/src/contrib/Archive/revss/revss_2.0.0.tar.gz"
  message("Installing revss v2.0.0 from CRAN archive into ", lib_path, " ...")
  install.packages(url, repos = NULL, type = "source", lib = lib_path, quiet = TRUE)
  if (!requireNamespace("revss", lib.loc = lib_path, quietly = TRUE)) {
    stop("Failed to install revss v2.0.0")
  }
  message("  OK")
}

# ---------------------------------------------------------------------------
# Capture function references from revss v2.0.0 (isolated load)
# ---------------------------------------------------------------------------
get_v2_fns <- function(lib_path) {
  ns <- loadNamespace("revss", lib.loc = lib_path)
  list(robLoc = ns$robLoc, robScale = ns$robScale, adm = ns$adm)
}

# ---------------------------------------------------------------------------
# Generic M-estimator benchmark (3 functions x sample sizes)
# ---------------------------------------------------------------------------
bench_m <- function(ns, loc_fn, scale_fn, adm_fn) {
  results <- vector("list", length(ns))
  for (i in seq_along(ns)) {
    n <- ns[i]
    set.seed(42L + n)
    x <- rnorm(n)
    iters <- get_iters(n)
    bm <- bench::mark(
      robLoc   = loc_fn(x),
      robScale = scale_fn(x),
      adm      = adm_fn(x),
      check = FALSE, iterations = iters
    )
    bm <- as.data.frame(bm)[, c("expression", "median")]
    bm$expr <- as.character(bm$expression)
    bm$n <- n
    results[[i]] <- bm[, c("expr", "n", "median")]
  }
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# Bonus: madn/admn benchmark (revss v3 vs robscale)
# ---------------------------------------------------------------------------
bench_bonus <- function(ns, madn_fn, admn_fn, mad_scaled_fn, adm_rs_fn) {
  results <- vector("list", length(ns))
  for (i in seq_along(ns)) {
    n <- ns[i]
    set.seed(42L + n)
    x <- rnorm(n)
    iters <- get_iters(n)
    bm <- bench::mark(
      `revss::madn`          = madn_fn(x),
      `revss::admn`          = admn_fn(x),
      `robscale::mad_scaled` = mad_scaled_fn(x),
      `robscale::adm`        = adm_rs_fn(x),
      check = FALSE, iterations = iters
    )
    bm <- as.data.frame(bm)[, c("expression", "median")]
    bm$expr <- as.character(bm$expression)
    bm$n <- n
    results[[i]] <- bm[, c("expr", "n", "median")]
  }
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# Build a three-way comparison table
# ---------------------------------------------------------------------------
build_comparison <- function(v2_df, v3_df, rs_df) {
  key <- function(df) paste(df$expr, df$n, sep = "|")
  v2_df$key <- key(v2_df)
  v3_df$key <- key(v3_df)
  rs_df$key <- key(rs_df)

  comp <- data.frame(
    expr     = v2_df$expr,
    n        = v2_df$n,
    t_v2     = as.numeric(v2_df$median),
    t_v3     = as.numeric(v3_df$median)[match(v2_df$key, v3_df$key)],
    t_rs     = as.numeric(rs_df$median)[match(v2_df$key, rs_df$key)],
    stringsAsFactors = FALSE
  )
  comp$v3_vs_v2 <- comp$t_v2 / comp$t_v3
  comp$rs_vs_v3 <- comp$t_v3 / comp$t_rs
  comp$rs_vs_v2 <- comp$t_v2 / comp$t_rs
  comp
}

# ---------------------------------------------------------------------------
# Build bonus comparison (revss v3 madn/admn vs robscale mad_scaled/adm)
# ---------------------------------------------------------------------------
build_bonus_comparison <- function(bonus_df) {
  revss_rows <- bonus_df[grepl("^revss::", bonus_df$expr), ]
  rs_rows    <- bonus_df[grepl("^robscale::", bonus_df$expr), ]

  revss_rows$pair <- ifelse(grepl("madn", revss_rows$expr), "mad", "adm")
  rs_rows$pair    <- ifelse(grepl("mad_scaled", rs_rows$expr), "mad", "adm")

  revss_rows$key <- paste(revss_rows$pair, revss_rows$n, sep = "|")
  rs_rows$key    <- paste(rs_rows$pair, rs_rows$n, sep = "|")

  idx <- match(revss_rows$key, rs_rows$key)
  data.frame(
    pair       = revss_rows$pair,
    n          = revss_rows$n,
    revss_fn   = revss_rows$expr,
    t_revss    = as.numeric(revss_rows$median),
    rs_fn      = rs_rows$expr[idx],
    t_rs       = as.numeric(rs_rows$median)[idx],
    stringsAsFactors = FALSE
  ) |>
    transform(rs_vs_revss = t_revss / t_rs)
}

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
fmt_time <- function(t_sec) {
  us <- t_sec * 1e6
  ifelse(us >= 1000, sprintf("%.1fms", us / 1000), sprintf("%.1fus", us))
}

print_main_table <- function(comp) {
  cat("\n", strrep("=", 90), "\n")
  cat("  revss v2.0.0 (R) vs revss v3.0.0 (Fortran) vs robscale (C++)\n")
  cat(strrep("=", 90), "\n\n")

  display <- data.frame(
    n       = comp$n,
    v2_time = fmt_time(comp$t_v2),
    v3_time = fmt_time(comp$t_v3),
    rs_time = fmt_time(comp$t_rs),
    `v3/v2` = sprintf("%.1fx", comp$v3_vs_v2),
    `rs/v3` = sprintf("%.1fx", comp$rs_vs_v3),
    `rs/v2` = sprintf("%.1fx", comp$rs_vs_v2),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  for (fn in unique(comp$expr)) {
    cat(sprintf("  --- %s ---\n", fn))
    rows <- comp$expr == fn
    print(display[rows, ], right = FALSE, row.names = FALSE)
    cat("\n")
  }
}

print_bonus_table <- function(bonus_comp) {
  cat("\n", strrep("=", 90), "\n")
  cat("  Bonus: revss v3.0.0 (madn/admn) vs robscale (mad_scaled/adm)\n")
  cat(strrep("=", 90), "\n\n")

  display <- data.frame(
    n          = bonus_comp$n,
    revss_fn   = bonus_comp$revss_fn,
    revss_time = fmt_time(bonus_comp$t_revss),
    rs_fn      = bonus_comp$rs_fn,
    rs_time    = fmt_time(bonus_comp$t_rs),
    speedup    = sprintf("%.1fx", bonus_comp$rs_vs_revss),
    stringsAsFactors = FALSE
  )

  for (p in unique(bonus_comp$pair)) {
    cat(sprintf("  --- %s ---\n", p))
    rows <- bonus_comp$pair == p
    print(display[rows, ], right = FALSE, row.names = FALSE)
    cat("\n")
  }
}

print_summary <- function(comp) {
  cat(strrep("=", 90), "\n")
  cat("  Summary (geometric mean speedup across all estimators & sizes)\n")
  cat(strrep("=", 90), "\n")
  cat(sprintf("  Fortran vs R   (v3/v2):  %.1fx\n", exp(mean(log(comp$v3_vs_v2)))))
  cat(sprintf("  C++ vs Fortran (rs/v3):  %.1fx\n", exp(mean(log(comp$rs_vs_v3)))))
  cat(sprintf("  C++ vs R       (rs/v2):  %.1fx\n", exp(mean(log(comp$rs_vs_v2)))))
  cat(strrep("=", 90), "\n")
}

# ===========================================================================
# Main
# ===========================================================================
main <- function() {
  start_time <- proc.time()
  message("Three-way benchmark: revss v2.0.0 vs v3.0.0 vs robscale\n")

  # 1. Install revss v2.0.0
  v2_lib <- tempfile("revss_v2_")
  install_revss_v2(v2_lib)
  v2_fns <- get_v2_fns(v2_lib)

  # 2. Get revss v3.0.0 functions (from normal library)
  library(revss)
  v3_fns <- list(
    robLoc   = revss::robLoc,
    robScale = revss::robScale,
    adm      = revss::adm
  )

  # 3. Get robscale functions
  library(robscale)
  rs_fns <- list(
    robLoc   = robscale::robLoc,
    robScale = robscale::robScale,
    adm      = robscale::adm
  )

  # 4. Run M-estimator benchmarks
  message("\nBenchmarking revss v2.0.0 (pure R) ...")
  v2_bm <- bench_m(sample_sizes, v2_fns$robLoc, v2_fns$robScale, v2_fns$adm)

  message("Benchmarking revss v3.0.0 (Fortran) ...")
  v3_bm <- bench_m(sample_sizes, v3_fns$robLoc, v3_fns$robScale, v3_fns$adm)

  message("Benchmarking robscale (C++) ...")
  rs_bm <- bench_m(sample_sizes, rs_fns$robLoc, rs_fns$robScale, rs_fns$adm)

  # 5. Bonus: madn/admn comparison
  message("Benchmarking bonus (madn/admn) ...")
  bonus_bm <- bench_bonus(
    sample_sizes,
    madn_fn       = revss::madn,
    admn_fn       = revss::admn,
    mad_scaled_fn = robscale::mad_scaled,
    adm_rs_fn     = robscale::adm
  )

  # 6. Build & print results
  comp <- build_comparison(v2_bm, v3_bm, rs_bm)
  print_main_table(comp)
  print_summary(comp)

  bonus_comp <- build_bonus_comparison(bonus_bm)
  print_bonus_table(bonus_comp)

  # 7. Save CSV
  elapsed <- (proc.time() - start_time)[["elapsed"]]
  cat(sprintf("\n  Total time: %.1f seconds\n\n", elapsed))

  out_file <- file.path("benchmarks",
                        paste0("bench_v2_v3_", format(Sys.Date(), "%Y%m%d"), ".csv"))
  write.csv(comp, out_file, row.names = FALSE)
  message("Results saved to ", out_file)
}

main()
