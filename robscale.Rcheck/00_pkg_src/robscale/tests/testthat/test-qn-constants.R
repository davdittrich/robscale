test_that("qn handles floating point variations in constant argument", {
  x <- rnorm(100)
  
  # Default call
  res_default <- qn(x)
  
  # Call with explicitly typed default constant
  res_explicit <- qn(x, constant = 2.2191)
  
  # Call with constant computed from simple arithmetic (which might drift depending on representation)
  val <- 1.10955 * 2 
  res_computed <- qn(x, constant = val)

  expect_equal(res_default, res_explicit)
  expect_equal(res_default, res_computed)
})
