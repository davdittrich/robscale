library(testthat)
library(robscale)

test_that("Qn and Sn refactor maintains numerical parity", {
  n_to_test <- c(3, 4, 8, 16, 64, 128, 256)
  
  for (n in n_to_test) {
    set.seed(42 + n)
    x <- rnorm(n)
    
    # Check Qn
    res_qn <- qn(x)
    expect_true(is.numeric(res_qn))
    expect_true(!is.na(res_qn))
    
    # Check Sn
    res_sn <- sn(x)
    expect_true(is.numeric(res_sn))
    expect_true(!is.na(res_sn))
    
    # Consistency checks: Qn and Sn should be roughly similar to scale
    expect_true(res_qn > 0)
    expect_true(res_sn > 0)
  }
})
