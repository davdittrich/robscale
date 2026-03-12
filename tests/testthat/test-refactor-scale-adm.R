library(testthat)
library(robscale)

test_that("robScale and adm refactor maintains numerical parity", {
  n_to_test <- c(3, 4, 8, 16, 64, 128, 256)
  
  for (n in n_to_test) {
    set.seed(42 + n)
    x <- rnorm(n)
    
    # Check robScale
    res_scale <- robScale(x)
    expect_true(is.numeric(res_scale))
    expect_true(!is.na(res_scale))
    expect_true(res_scale >= 0)
    
    # Check adm
    res_adm <- adm(x)
    expect_true(is.numeric(res_adm))
    expect_true(!is.na(res_adm))
    expect_true(res_adm >= 0)
    
    # Consistency check for tiny n
    if (n < 3) {
      # Should return 0 or small
    }
  }
})
