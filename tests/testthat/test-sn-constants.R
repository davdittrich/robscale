test_that("sn handles floating point variations in constant argument", {
  x <- rnorm(100)

  # Default call
  res_default <- sn(x)

  # Call with explicitly typed full-precision default constant
  res_explicit <- sn(x, constant = 1.19259855312321)

  # Call with constant computed from simple arithmetic
  val <- 0.596299276561605 * 2
  res_computed <- sn(x, constant = val)

  expect_equal(res_default, res_explicit)
  expect_equal(res_default, res_computed)
})
