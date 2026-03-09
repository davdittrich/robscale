test_that("robScale fallback = 'adm' (default) handles MAD collapse", {
  # Degenerate input: MAD will be 0
  x <- c(1, 1, 1, 1, 5)
  
  # Default behavior: fallback to adm()
  expect_no_error(s <- robScale(x))
  expect_true(s > 0)
  expect_equal(s, adm(x))
  
  # Explicit fallback = "adm"
  expect_equal(robScale(x, fallback = "adm"), adm(x))
})

test_that("robScale fallback = 'na' returns NA on MAD collapse", {
  x <- c(1, 1, 1, 1, 5)
  
  # fallback = "na" should return NA
  expect_true(is.na(robScale(x, fallback = "na")))
})

test_that("robScale fallback = 'na' returns NA for very small samples if MAD collapses", {
  x <- c(1, 1) # n=2 < minobs=4
  expect_true(is.na(robScale(x, fallback = "na")))
  
  # n=2, values different -> MAD > 0
  x2 <- c(1, 2)
  expect_false(is.na(robScale(x2, fallback = "na")))
})

test_that("robScale handles implbound correctly", {
  # MAD is very small but > 0
  x <- c(1, 1, 1, 1, 1.000001)
  
  # With implbound > 0, it should fallback
  expect_equal(robScale(x, implbound = 1e-4), adm(x))
  
  # With implbound = 0, it should not fallback (unless truly 0)
  # But it might still struggle with literal 0 if it's below something else.
  # Here MAD is ~1e-6.
  expect_true(robScale(x, implbound = 0) > 0)
})
