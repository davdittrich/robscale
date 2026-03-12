test_that("robscale:::get_qnsn_config() returns cache-derived thresholds", {
  config <- robscale:::get_qnsn_config()

  # Grain size should be >= 512 and <= 8192
  expect_gte(config$grain_size, 512L)
  expect_lte(config$grain_size, 8192L)

  # Parallel thresholds should be >= 4096
  expect_gte(config$qn_parallel_threshold, 4096L)
  expect_gte(config$sn_parallel_threshold, 4096L)

  # Sort TBB threshold should be >= 4096
  expect_gte(config$sort_tbb_threshold, 4096L)

  # Sort boost threshold should be reasonable
  expect_gte(config$sort_boost_threshold, 64L)
  expect_lte(config$sort_boost_threshold, 8192L)

  # Fixed thresholds should remain constant
  expect_equal(config$qn_exact_threshold, 64L)
})

test_that("generated Makevars has no dead macros", {
  skip_if(!file.exists("src/Makevars"), "Not in source tree")
  content <- paste(readLines("src/Makevars"), collapse = " ")
  expect_false(grepl("ROBSCALE_FAST_MODE", content))
  expect_false(grepl("FASTQNSN_L2_CACHE_BYTES", content))
  expect_false(grepl("FASTQNSN_CACHE_LINE_BYTES", content))
})

test_that("SIMD level matches platform", {
  config <- robscale:::get_qnsn_config()
  arch <- Sys.info()[["machine"]]
  if (arch %in% c("arm64", "aarch64")) {
    # ARM64: NEON is always available (simd_level = 4 = Neon)
    expect_equal(config$simd_level, 4L)
  } else if (arch %in% c("x86_64", "amd64")) {
    # x86_64: at minimum SSE4.2 (simd_level >= 1)
    expect_gte(config$simd_level, 1L)
  }
})
