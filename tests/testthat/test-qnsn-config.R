## tests/testthat/test-qnsn-config.R
## Tests for WU-PORT-1 portability fixes in qnsn_hardware_info.h

test_that("get_qnsn_config() returns a list", {
  cfg <- robscale:::get_qnsn_config()
  expect_type(cfg, "list")
})

test_that("get_qnsn_config() numeric fields are finite, positive, and < 2^31", {
  cfg <- robscale:::get_qnsn_config()
  numeric_fields <- c(
    "l2_cache_size", "l2_per_core", "num_logical_cores", "num_physical_cores",
    "qn_parallel_threshold", "sn_parallel_threshold", "sort_tbb_threshold",
    "grain_size"
  )
  for (fld in numeric_fields) {
    val <- cfg[[fld]]
    expect_true(is.finite(val), info = paste(fld, "is not finite"))
    expect_true(val > 0,        info = paste(fld, "is not > 0"))
    expect_true(val < 2^31,     info = paste(fld, "is not < 2^31"))
  }
})

test_that("get_qnsn_config()$num_logical_cores >= 1L", {
  cfg <- robscale:::get_qnsn_config()
  expect_gte(cfg$num_logical_cores, 1L)
})
