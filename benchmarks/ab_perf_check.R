#!/usr/bin/env Rscript --vanilla
# A/B performance comparison: old (installed) vs new (source) robscale
# Installs the current source into a temp library, then benchmarks both.

library(bench)

cat("=== A/B Performance Test ===\n")
cat("Platform:", R.version$platform, "\n")
cat("R version:", R.version.string, "\n\n")

# --- Build new version into temp library ---
source_path <- if (file.exists("DESCRIPTION")) "." else ".."
new_lib <- tempfile("robscale_new_")
dir.create(new_lib)

# Copy source to clean temp dir (avoid stale .o files)
tmp_src <- tempfile("src_")
dir.create(tmp_src)
for (f in c("DESCRIPTION", "NAMESPACE", "configure", "cleanup")) {
  file.copy(file.path(source_path, f), tmp_src)
}
for (d in c("R", "src", "inst")) {
  file.copy(file.path(source_path, d), tmp_src, recursive = TRUE)
}
# Purge compiled objects
unlink(list.files(file.path(tmp_src, "src"),
                  pattern = "\\.o$|\\.so$|\\.dll$",
                  full.names = TRUE))

cat("Installing new version from source...\n")
install_result <- system2(
  R.home("bin/R"),
  args = c("CMD", "INSTALL", tmp_src,
           paste0("--library=", new_lib),
           "--no-multiarch", "--no-test-load"),
  stdout = TRUE, stderr = TRUE
)
unlink(tmp_src, recursive = TRUE)
if (!is.null(attr(install_result, "status"))) {
  cat(paste(install_result, collapse = "\n"), "\n")
  stop("Installation of new version failed")
}
cat("New version installed to:", new_lib, "\n\n")

# --- Benchmark function ---
run_bench <- function(lib_path, label) {
  callr::r(function(lib_path, label) {
    if (!is.null(lib_path)) {
      .libPaths(c(lib_path, .libPaths()))
    }
    library(robscale)
    library(bench)

    cat(sprintf("--- %s (robscale %s from %s) ---\n",
                label, packageVersion("robscale"),
                dirname(system.file("DESCRIPTION", package = "robscale"))))

    n_grid <- c(16, 64, 256, 1000, 4096, 10000, 50000, 100000, 500000, 1000000)

    get_iters <- function(n) {
      if (n <= 128) 5000L
      else if (n <= 2048) 1000L
      else if (n <= 16384) 200L
      else if (n <= 100000) 50L
      else 10L
    }

    iqr_res <- bench::press(n = n_grid, {
      set.seed(42 + n); x <- rnorm(n)
      bench::mark(iqr = robscale::iqr_scaled(x),
                  check = FALSE, min_iterations = get_iters(n), min_time = 0.5)
    })

    mad_res <- bench::press(n = n_grid, {
      set.seed(42 + n); x <- rnorm(n)
      bench::mark(mad = robscale::mad_scaled(x),
                  check = FALSE, min_iterations = get_iters(n), min_time = 0.5)
    })

    gmd_res <- bench::press(n = n_grid, {
      set.seed(42 + n); x <- rnorm(n)
      bench::mark(gmd = robscale::gmd(x),
                  check = FALSE, min_iterations = get_iters(n), min_time = 0.5)
    })

    adm_res <- bench::press(n = n_grid, {
      set.seed(42 + n); x <- rnorm(n)
      bench::mark(adm = robscale::adm(x),
                  check = FALSE, min_iterations = get_iters(n), min_time = 0.5)
    })

    ensemble_res <- bench::press(n = c(32, 128, 512), {
      set.seed(42 + n); x <- rnorm(n)
      bench::mark(ensemble = robscale::scale_robust(x, method = "ensemble"),
                  check = FALSE, min_iterations = 20L, min_time = 1.0)
    })

    list(
      iqr = iqr_res[, c("n", "median", "mem_alloc")],
      mad = mad_res[, c("n", "median", "mem_alloc")],
      gmd = gmd_res[, c("n", "median", "mem_alloc")],
      adm = adm_res[, c("n", "median", "mem_alloc")],
      ensemble = ensemble_res[, c("n", "median", "mem_alloc")]
    )
  }, args = list(lib_path = lib_path, label = label), show = TRUE)
}

# --- Run OLD version (system library) ---
cat("Running OLD (installed) version...\n")
old <- run_bench(NULL, "OLD")

# --- Run NEW version (temp library) ---
cat("\nRunning NEW (optimized) version...\n")
new <- run_bench(new_lib, "NEW")

# --- Compare ---
compare <- function(old_df, new_df, name) {
  merged <- merge(old_df, new_df, by = "n", suffixes = c("_old", "_new"))
  merged$speedup <- as.numeric(merged$median_old) / as.numeric(merged$median_new)
  cat(sprintf("\n=== %s Speedup (old/new) ===\n", name))
  for (i in seq_len(nrow(merged))) {
    cat(sprintf("  n=%7d  old=%10s  new=%10s  speedup=%.2fx  mem_old=%s  mem_new=%s\n",
                merged$n[i],
                format(merged$median_old[i]),
                format(merged$median_new[i]),
                merged$speedup[i],
                format(merged$mem_alloc_old[i]),
                format(merged$mem_alloc_new[i])))
  }
  cat(sprintf("  MEDIAN SPEEDUP: %.2fx\n", median(merged$speedup)))
}

compare(old$iqr, new$iqr, "IQR (P1: incremental offset)")
compare(old$mad, new$mad, "MAD (P4: fused median-then-MAD)")
compare(old$gmd, new$gmd, "GMD (P2: SIMD hints)")
compare(old$adm, new$adm, "ADM (P2: SIMD hints)")
compare(old$ensemble, new$ensemble, "Ensemble (compounded)")

# Cleanup
unlink(new_lib, recursive = TRUE)
cat("\nDone.\n")
