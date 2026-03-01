tol <- sqrt(.Machine$double.eps)

if (requireNamespace("revss", quietly = TRUE)) {
  set.seed(42)
  for (n in 3:20) {
    for (rep in 1:100) {
      x <- runif(n, -100, 100)

      # ADM
      expect_equal(
        robscale::adm(x), revss::adm(x),
        tolerance = tol,
        info = paste0("adm n=", n, " rep=", rep)
      )

      # robLoc
      expect_equal(
        robscale::robLoc(x), revss::robLoc(x),
        tolerance = tol,
        info = paste0("robLoc n=", n, " rep=", rep)
      )

      # robScale
      expect_equal(
        robscale::robScale(x), revss::robScale(x),
        tolerance = tol,
        info = paste0("robScale n=", n, " rep=", rep)
      )
    }
  }
}
