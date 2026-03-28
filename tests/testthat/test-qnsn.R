test_that("sn() works for basic cases", {
  x <- c(1, 2, 4, 8, 16)
  expect_no_error(s <- sn(x))
  expect_true(is.numeric(s))
  expect_true(s > 0)
  
  # Known value for n=5
  # sorted x: 1, 2, 4, 8, 16
  # |x_i - x_j| for each i:
  # i=1 (1): |1-1|=0, |1-2|=1, |1-4|=3, |1-8|=7, |1-16|=15. lowmed: 1
  # i=2 (2): |2-1|=1, |2-2|=0, |2-4|=2, |2-8|=6, |2-16|=14. lowmed: 2
  # i=3 (4): |4-1|=3, |4-2|=2, |4-4|=0, |4-8|=4, |4-16|=12. lowmed: 3
  # i=4 (8): |8-1|=7, |8-2|=6, |8-4|=4, |8-8|=0, |8-16|=8.  lowmed: 6
  # i=5 (16):|16-1|=15,|16-2|=14,|16-4|=12,|16-8|=8,|16-16|=0.lowmed: 12
  # inner medians: 1, 2, 3, 6, 12. lowmed of inner: 3
  # result: 3 * 1.1926 * factor_sn(5=1.34857) = 3 * 1.6083... = 4.825
  expect_equal(sn(x, finite.corr = FALSE), 3 * 1.19259855312321)
})

test_that("qn() works for basic cases", {
  x <- c(1, 2, 4, 8, 16)
  expect_no_error(q <- qn(x))
  expect_true(is.numeric(q))
  expect_true(q > 0)
  
  # pairwise diffs: 1, 3, 7, 15, 2, 6, 14, 4, 12, 8
  # sorted: 1, 2, 3, 4, 6, 7, 8, 12, 14, 15
  # n=5, h=3, k = 3*(3-1)/2 = 3
  # 3rd element: 3
  # result: 3 * 2.2191 * factor_qn(5=0.84412)
  expect_equal(qn(x, finite.corr = FALSE), 3 * 2.21914446598508)
})

test_that("sn() and qn() handle NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  # R7: now errors on NA with na.rm=FALSE (consistent with all other wrappers)
  expect_error(sn(x), "NAs")
  expect_error(qn(x), "NAs")

  expect_equal(sn(x, na.rm = TRUE), sn(x[!is.na(x)]))
  expect_equal(qn(x, na.rm = TRUE), qn(x[!is.na(x)]))
})

test_that("sn() and qn() handle small samples", {
  expect_true(is.na(sn(1)))
  expect_true(is.na(qn(1)))
  
  expect_true(is.numeric(sn(c(1, 2))))
  expect_true(is.numeric(qn(c(1, 2))))
})

test_that("sn() and qn() handle integer inputs", {
  x <- c(1L, 2L, 4L, 8L, 16L)
  expect_equal(sn(x, finite.corr = FALSE), 3 * 1.19259855312321)
  expect_equal(qn(x, finite.corr = FALSE), 3 * 2.21914446598508)
})
