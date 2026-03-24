test_that("mad micro-buffer: n in 65:128 uses stack path (correct output)", {
  for (n in c(65L, 80L, 100L, 128L)) {
    set.seed(42)
    x <- rnorm(n)
    # mad_scaled(constant=1) should match stats::mad(x, constant=1) to 8 ULP
    expect_equal(
      mad_scaled(x, constant = 1),
      stats::mad(x, constant = 1),
      tolerance = 8 * .Machine$double.eps,
      label = paste0("n=", n)
    )
  }
})

test_that("mad micro-buffer: n in 65:128 with user-supplied center matches stats::mad", {
  for (n in c(65L, 80L, 100L, 128L)) {
    set.seed(99)
    x <- rnorm(n)
    center_val <- median(x)
    # mad_scaled with center supplied should match stats::mad(x, center=center_val, constant=1)
    expect_equal(
      mad_scaled(x, center = center_val, constant = 1),
      stats::mad(x, center = center_val, constant = 1),
      tolerance = 8 * .Machine$double.eps,
      label = paste0("n=", n, " (center path)")
    )
  }
})
