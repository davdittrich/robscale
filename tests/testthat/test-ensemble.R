test_that("ensemble is deterministic", {
  set.seed(42)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 200L)
  e2 <- robscale:::cpp_scale_ensemble(x, 200L)
  expect_identical(e1, e2)
})

test_that("ensemble estimates sigma under normality", {
  set.seed(42)
  y <- rnorm(5000)
  expect_equal(robscale:::cpp_scale_ensemble(y, 200L), 1.0, tolerance = 0.05)
})

test_that("ensemble is more robust than sd for contaminated data", {
  set.seed(42)
  x_clean <- rnorm(20)
  x_dirty <- c(x_clean, 100)  # one extreme outlier
  ens_clean <- robscale:::cpp_scale_ensemble(x_clean, 200L)
  ens_dirty <- robscale:::cpp_scale_ensemble(x_dirty, 200L)
  sd_dirty  <- sd(x_dirty)

  # Ensemble should be much less affected than sd
  expect_true(abs(ens_dirty - ens_clean) < abs(sd_dirty - sd(x_clean)))
})

test_that("ensemble handles edge cases", {
  expect_true(is.na(robscale:::cpp_scale_ensemble(1, 200L)))
  expect_true(is.numeric(robscale:::cpp_scale_ensemble(c(1, 2), 200L)))
  expect_true(robscale:::cpp_scale_ensemble(c(1, 2), 200L) > 0)
})

test_that("ensemble returns finite for small samples", {
  for (n in 2:10) {
    set.seed(n)
    x <- rnorm(n)
    ens <- robscale:::cpp_scale_ensemble(x, 200L)
    expect_true(is.finite(ens), label = paste("n =", n))
    expect_true(ens > 0, label = paste("n =", n))
  }
})

test_that("ensemble handles data with ties", {
  # Most robust estimators collapse, but ensemble should still return a value
  val <- robscale:::cpp_scale_ensemble(c(5, 5, 5, 5, 6), 200L)
  expect_true(is.finite(val))
})

# --- Snapshot and determinism stress tests ---

test_that("ensemble snapshot values are stable", {
  # Tolerance 1e-12: accommodates ULP-level differences across platforms
  # (ARM64 vs x86_64) and optimization levels (-O0 vs -O3 with SIMD).
  # Determinism tests below use expect_identical for parallelism safety.
  set.seed(1)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(5), 200L),
    1.033903437859763796, tolerance = 1e-12
  )
  set.seed(42)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(10), 200L),
    0.8309512014065341123, tolerance = 1e-12
  )
  set.seed(123)
  expect_equal(
    robscale:::cpp_scale_ensemble(rnorm(500), 200L),
    0.9615157792473906229, tolerance = 1e-12
  )
})

test_that("ensemble is deterministic over 50 repeats (n=10)", {
  set.seed(99)
  x <- rnorm(10)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(50)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble is deterministic for various n", {
  for (n in c(2, 5, 10, 15, 50, 200)) {
    set.seed(n)
    x <- rnorm(n)
    e1 <- robscale:::cpp_scale_ensemble(x, 200L)
    e2 <- robscale:::cpp_scale_ensemble(x, 200L)
    expect_identical(e1, e2, label = paste("n =", n))
  }
})

test_that("ensemble is deterministic with large n_boot", {
  set.seed(7)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 500L)
  e2 <- robscale:::cpp_scale_ensemble(x, 500L)
  expect_identical(e1, e2)
})

test_that("ensemble is deterministic with large n", {
  set.seed(11)
  x <- rnorm(500)
  e1 <- robscale:::cpp_scale_ensemble(x, 200L)
  e2 <- robscale:::cpp_scale_ensemble(x, 200L)
  expect_identical(e1, e2)
})

test_that("ensemble determinism at threshold boundary (n*n_boot=10000)", {
  set.seed(33)
  x <- rnorm(50)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(20)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble determinism above threshold (n=100)", {
  set.seed(44)
  x <- rnorm(100)
  ref <- robscale:::cpp_scale_ensemble(x, 200L)
  for (i in seq_len(20)) {
    expect_identical(robscale:::cpp_scale_ensemble(x, 200L), ref)
  }
})

test_that("ensemble determinism high n_boot at threshold (n=10, n_boot=1000)", {
  set.seed(55)
  x <- rnorm(10)
  e1 <- robscale:::cpp_scale_ensemble(x, 1000L)
  e2 <- robscale:::cpp_scale_ensemble(x, 1000L)
  expect_identical(e1, e2)
})

# --- Ensemble CI tests ---

test_that("scale_robust(ci=TRUE) returns robscale_ensemble_ci class", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = TRUE)
  expect_s3_class(res, "robscale_ensemble_ci")
  expect_true(is.numeric(res$estimate))
  expect_named(res$ci, c("lower", "upper"))
  expect_equal(res$level, 0.95)
  expect_equal(res$method, "ensemble")
})

test_that("ensemble estimate is within its bootstrap CI", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = TRUE)
  expect_true(res$ci[["lower"]] <= res$estimate)
  expect_true(res$ci[["upper"]] >= res$estimate)
})

test_that("all 7 per-estimator rows present with both CI types", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = TRUE)
  df <- res$estimators
  expect_equal(nrow(df), 7L)
  expected_names <- c("sd_c4", "gmd", "mad_scaled", "iqr_scaled",
                      "sn", "qn", "robScale")
  expect_equal(df$estimator, expected_names)
  expect_true(all(is.finite(df$analytical_lower)))
  expect_true(all(is.finite(df$analytical_upper)))
  expect_true(all(is.finite(df$boot_lower)))
  expect_true(all(is.finite(df$boot_upper)))
})

test_that("ensemble weights sum to 1", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = TRUE)
  expect_equal(sum(res$estimators$weight), 1.0, tolerance = 1e-10)
})

test_that("ensemble ci: level=0.99 wider than level=0.90", {
  set.seed(42)
  x <- rnorm(10)
  r90 <- scale_robust(x, ci = TRUE, level = 0.90)
  r99 <- scale_robust(x, ci = TRUE, level = 0.99)
  w90 <- r90$ci[["upper"]] - r90$ci[["lower"]]
  w99 <- r99$ci[["upper"]] - r99$ci[["lower"]]
  expect_true(w99 > w90)
})

test_that("ensemble ci is deterministic", {
  set.seed(42)
  x <- rnorm(10)
  r1 <- scale_robust(x, ci = TRUE)
  r2 <- scale_robust(x, ci = TRUE)
  expect_identical(r1$estimate, r2$estimate)
  expect_identical(r1$ci, r2$ci)
})

test_that("ci=FALSE still returns scalar from scale_robust", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = FALSE)
  expect_true(is.numeric(res))
  expect_length(res, 1)
})

test_that("ci=TRUE with single method returns robscale_ci", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, method = "qn", ci = TRUE)
  expect_s3_class(res, "robscale_ci")
  expect_equal(res$method, "qn")
})

test_that("ci=TRUE with auto_switch returns robscale_ci for gmd", {
  set.seed(42)
  x <- rnorm(50)  # n >= 20 triggers auto_switch
  res <- scale_robust(x, ci = TRUE)
  expect_s3_class(res, "robscale_ci")
  expect_equal(res$method, "gmd")
})

test_that("BCa z0 and acc values are finite and reasonable", {
  set.seed(42)
  x <- rnorm(10)
  raw <- robscale:::cpp_scale_ensemble_ci(x, 200L, 0.95, 0L)
  expect_true(all(is.finite(raw$z0)))
  expect_true(all(is.finite(raw$acc)))
  # z0 should be close to 0 for symmetric data
  expect_true(all(abs(raw$z0) < 3))
  # acc should be small
  expect_true(all(abs(raw$acc) < 1))
})

test_that("boot_method tiers work correctly", {
  set.seed(42)
  x <- rnorm(10)
  bca <- scale_robust(x, ci = TRUE, boot_method = "bca")
  pct <- scale_robust(x, ci = TRUE, boot_method = "percentile")
  par <- scale_robust(x, ci = TRUE, boot_method = "parametric")
  expect_equal(bca$boot_method, "bca")
  expect_equal(pct$boot_method, "percentile")
  expect_equal(par$boot_method, "parametric")
  # All should produce valid CIs
  expect_true(bca$ci[["lower"]] < bca$ci[["upper"]])
  expect_true(pct$ci[["lower"]] < pct$ci[["upper"]])
  expect_true(par$ci[["lower"]] < par$ci[["upper"]])
})

test_that("print.robscale_ensemble_ci produces output", {
  set.seed(42)
  x <- rnorm(10)
  res <- scale_robust(x, ci = TRUE)
  expect_output(print(res), "Ensemble estimate:")
  expect_output(print(res), "bootstrap CI")
})
