test_that("get_qnsn_config exposes hardware diagnostics", {
  cfg <- robscale:::get_qnsn_config()
  for (f in c("l2_per_core", "num_physical_cores", "has_tuned_sort_thresholds"))
    expect_true(f %in% names(cfg), info = paste("Missing:", f))
})

test_that("l2_per_core is between 128KB and 16MB", {
  cfg <- robscale:::get_qnsn_config()
  expect_gte(cfg$l2_per_core, 131072L)    # 128KB min
  expect_lte(cfg$l2_per_core, 16777216L)  # 16MB max
})

test_that("sort_tbb_threshold >= 8192 on any modern CPU", {
  cfg <- robscale:::get_qnsn_config()
  expect_gte(cfg$sort_tbb_threshold, 8192L)
})

test_that("sort thresholds are hardware-derived, not tuner-overridden", {
  cfg <- robscale:::get_qnsn_config()
  expect_false(cfg$has_tuned_sort_thresholds)
})

test_that("thresholds are derived from l2_per_core", {
  cfg <- robscale:::get_qnsn_config()
  l2 <- cfg$l2_per_core
  expect_equal(cfg$sort_tbb_threshold, max(4096L, as.integer(l2 / 16)))
  expect_equal(cfg$qn_parallel_threshold, max(4096L, as.integer(l2 / 16)))
  expect_equal(cfg$sn_parallel_threshold, max(4096L, as.integer(l2 / 16)))
  expect_equal(cfg$grain_size, min(8192L, max(512L, as.integer(l2 / 32))))
})

test_that("Qn/Sn correct across threshold boundaries", {
  set.seed(42)
  for (n in c(10, 64, 512, 4096, 8192, 16384)) {
    x <- rnorm(n)
    expect_true(is.finite(qn(x)) && qn(x) > 0)
    expect_true(is.finite(sn(x)) && sn(x) > 0)
  }
})
