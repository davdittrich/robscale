test_that("robscale matches revss baseline", {
  skip_if_not_installed("revss")
  tol <- sqrt(.Machine$double.eps)
  set.seed(42)
  for (n in 3:20) { 
    for (rep in 1:100) {
      x <- runif(n, -100, 100)
      expect_equal(robscale::adm(x), revss::adm(x), tolerance = tol)
      expect_equal(robscale::robLoc(x), revss::robLoc(x), tolerance = tol)
      expect_equal(robscale::robScale(x), revss::robScale(x), tolerance = tol)
    }
  }
})
