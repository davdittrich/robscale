test_that("qn handles floating point variations in constant argument", {
  x <- rnorm(100)

  # Default call
  res_default <- qn(x)

  # Call with explicitly typed full-precision default constant
  res_explicit <- qn(x, constant = 2.21914446598508)

  # Call with constant computed from simple arithmetic
  val <- 1.10957223299254 * 2
  res_computed <- qn(x, constant = val)

  expect_equal(res_default, res_explicit)
  expect_equal(res_default, res_computed)
})
