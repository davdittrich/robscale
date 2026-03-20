test_that("get_sn_factor handles n < 2 safely", {
  expect_equal(robscale:::C_get_sn_factor(0), 0.0)
  expect_equal(robscale:::C_get_sn_factor(1), 0.0)
  expect_gt(robscale:::C_get_sn_factor(2), 0.0)
})

test_that("get_qn_factor handles n < 2 safely", {
  expect_equal(robscale:::C_get_qn_factor(0), 0.0)
  expect_equal(robscale:::C_get_qn_factor(1), 0.0)
  expect_gt(robscale:::C_get_qn_factor(2), 0.0)
})
