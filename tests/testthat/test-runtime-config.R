test_that("get_qnsn_config exposes hardware diagnostics", {
  cfg <- robscale:::get_qnsn_config()
  for (f in c("l2_per_core", "num_physical_cores", "has_tuned_sort_thresholds"))
    expect_true(f %in% names(cfg), info = paste("Missing:", f))
})

test_that("l2_per_core is between 64KB and 16MB", {
  cfg <- robscale:::get_qnsn_config()
  expect_gte(cfg$l2_per_core, 65536L)     # 64KB min (sanity floor)
  expect_lte(cfg$l2_per_core, 16777216L)  # 16MB max
})

test_that("sort_tbb_threshold >= 4096 on any modern CPU", {
  cfg <- robscale:::get_qnsn_config()
  expect_gte(cfg$sort_tbb_threshold, 4096L)
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

test_that("pdq_robscale_threshold present and in valid range", {
  cfg <- robscale:::get_qnsn_config()
  expect_true("pdq_robscale_threshold" %in% names(cfg))
  thr <- cfg$pdq_robscale_threshold
  expect_gte(thr, 2048L)
  expect_lte(thr, 500000L)
  # Verify formula: max(2048, L2 / (sizeof(double) * 2)) == max(2048, L2 / 16)
  expect_equal(thr, max(2048L, as.integer(cfg$l2_per_core / 16)))
})

test_that("pdq_lowmedian_threshold present and in valid range", {
  cfg <- robscale:::get_qnsn_config()
  expect_true("pdq_lowmedian_threshold" %in% names(cfg))
  thr <- cfg$pdq_lowmedian_threshold
  expect_gte(thr, 2048L)
  expect_lte(thr, 500000L)
  # Verify formula: max(2048, L2 / (sizeof(double) * 2)) == max(2048, L2 / 16)
  expect_equal(thr, max(2048L, as.integer(cfg$l2_per_core / 16)))
})

test_that("pdq_qn_final_threshold present and in valid range", {
  cfg <- robscale:::get_qnsn_config()
  expect_true("pdq_qn_final_threshold" %in% names(cfg))
  thr <- cfg$pdq_qn_final_threshold
  expect_gte(thr, 2048L)
  expect_lte(thr, 500000L)
  # Verify formula: max(2048, L2 / (sizeof(double) * 4)) == max(2048, L2 / 32)
  expect_equal(thr, max(2048L, as.integer(cfg$l2_per_core / 32)))
})
