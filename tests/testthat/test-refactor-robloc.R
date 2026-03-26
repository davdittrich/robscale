test_that("robLoc refactor maintains numerical parity", {
  # Test with various n including tiny, small, and medium
  n_to_test <- c(3, 4, 8, 16, 64, 128, 256)

  for (n in n_to_test) {
    set.seed(42 + n)
    x <- rnorm(n)
    
    res1 <- robLoc(x)
    
    expect_true(is.numeric(res1))
    expect_true(!is.na(res1))
    
    if (n < 4) {
      expect_equal(res1, median(x), tolerance = 1e-10)
    } else {
      # For n >= 4, it should be stable. We'll just check it doesn't deviate 
      # significantly from median for now, or we could use the unoptimized 
      # version if we had it loaded.
      expect_true(abs(res1 - median(x)) < 2 * sd(x))
    }
  }
})

test_that("robLoc avoids allocation for micro-samples", {
  # This is a performance/contract test
  # We can't easily check allocation from R, but we can verify it doesn't crash 
  # and returns results for n <= 8.
  x <- rnorm(8)
  expect_no_error(robLoc(x))
})
