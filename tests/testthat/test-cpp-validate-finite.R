# WU-VAL-1: C++ validate_finite at all entry points
# These tests call internal C++ exports DIRECTLY (bypassing R wrappers).
# Before the fix, non-finite inputs return garbage/NA instead of erroring.

test_that("C++ rob_scale_impl errors on Inf", {
  expect_error(robscale:::C_rob_scale_fast(c(1, Inf, 3, 4, 5)), "finite")
})

test_that("C++ rob_loc_impl errors on Inf", {
  expect_error(robscale:::rob_loc_impl(c(1, Inf, 3, 4, 5),
    FALSE, 0.0, 80L, sqrt(.Machine$double.eps)), "finite")
})

test_that("C++ sd_c4_impl errors on Inf", {
  expect_error(robscale:::sd_c4_impl(c(1, Inf, 3, 4, 5)), "finite")
})

test_that("C++ adm_impl_auto errors on Inf", {
  expect_error(robscale:::adm_impl_auto(c(1, Inf, 3, 4, 5),
    1.2533141373155001), "finite")
})

test_that("C++ gmd_impl errors on Inf", {
  expect_error(robscale:::gmd_impl(c(1, Inf, 3, 4, 5),
    0.886226925452758), "finite")
})

test_that("C++ mad_impl_auto errors on Inf", {
  expect_error(robscale:::mad_impl_auto(c(1, Inf, 3, 4, 5),
    1.4826022185056), "finite")
})

test_that("C++ iqr_impl errors on Inf", {
  expect_error(robscale:::iqr_impl(c(1, Inf, 3, 4, 5),
    0.7413011092528), "finite")
})

test_that("C++ C_qn_fast errors on Inf (was R_NaReal)", {
  expect_error(robscale:::C_qn_fast(c(1, Inf, 3, 4, 5)), "finite")
})

test_that("C++ C_sn_fast errors on Inf (was R_NaReal)", {
  expect_error(robscale:::C_sn_fast(c(1, Inf, 3, 4, 5)), "finite")
})

test_that("C++ sd_c4_impl errors on NA", {
  expect_error(robscale:::sd_c4_impl(c(1, NA, 3, 4, 5)), "NAs")
})

test_that("C++ cpp_scale_ensemble errors on Inf", {
  expect_error(robscale:::cpp_scale_ensemble(c(1, Inf, 3, 4, 5), 100L), "finite")
})

test_that("C++ robScale with loc errors on Inf", {
  # has_loc=TRUE path in rob_scale_core
  expect_error(robScale(c(1, Inf, 3, 4, 5), loc = 0), "finite")
})

test_that("mixed na.rm + Inf: NA stripped then Inf caught by C++", {
  expect_error(qn(c(1, NA, Inf, 3, 4), na.rm = TRUE), "finite")
})
