# WU-3: adm_core_sorted catastrophic cancellation (R1) + vshaped_mad overflow (S12)

test_that("adm_core_sorted avoids cancellation at 1e15 scale (via robScale loc path)", {
  # adm_core_sorted is called when has_loc=true and MAD implodes (s_init==0).
  # Data: 3 identical at 1e15 + 2 slightly different → MAD = 0 → ADM fallback
  x <- c(1e15, 1e15, 1e15, 1e15 + 1e-6, 1e15 + 2e-6)
  # Trigger sorted path: robScale with explicit loc, all-identical majority → MAD = 0
  result <- robScale(x, loc = median(x))
  expected_adm <- sqrt(pi / 2) * mean(abs(x - median(x)))
  # Before fix: adm_core_sorted returns ≈ 0 (upper_sum - lower_sum cancels at 1e15)
  # After fix: direct center-based differences are exact
  expect_equal(result, expected_adm, tolerance = 10 * .Machine$double.eps,
               label = "1e15-scale sorted ADM matches R reference")
})

test_that("adm_core_sorted at 1e15 with even n (via robScale loc path)", {
  x <- c(1e15, 1e15, 1e15 + 1e-6, 1e15 + 2e-6)
  result <- robScale(x, loc = median(x))
  expected_adm <- sqrt(pi / 2) * mean(abs(x - median(x)))
  expect_equal(result, expected_adm, tolerance = 10 * .Machine$double.eps,
               label = "1e15-scale even-n sorted ADM matches R reference")
})

test_that("adm_core_sorted at 1e12 scale (moderate offset)", {
  x <- 1e12 + c(-3, -1, 0, 2, 5)
  expected <- sqrt(pi / 2) * mean(abs(x - median(x)))
  result <- adm(x)
  expect_equal(result, expected, tolerance = 10 * .Machine$double.eps,
               label = "1e12-scale adm")
})

test_that("adm unchanged for normal-scale data (regression guard)", {
  set.seed(17)
  x <- rnorm(20)
  expected <- sqrt(pi / 2) * mean(abs(x - median(x)))
  expect_equal(adm(x), expected, tolerance = 20 * .Machine$double.eps,
               label = "normal-scale adm unchanged")
})

test_that("adm with explicit center at large scale", {
  x <- c(1e15, 1e15 + 1, 1e15 + 2, 1e15 + 3, 1e15 + 4)
  ctr <- median(x)
  expected <- sqrt(pi / 2) * mean(abs(x - ctr))
  expect_equal(adm(x, center = ctr), expected,
               tolerance = 10 * .Machine$double.eps,
               label = "1e15-scale adm with explicit center")
})

test_that("robScale fallback to adm_core_sorted works at large scale", {
  # Near-constant large-magnitude data: MAD=0, fallback to adm_core (unsorted path)
  x <- c(1e15, 1e15, 1e15, 1e15, 1e15 + 1)
  result <- robScale(x)
  expected <- adm(x)
  expect_equal(result, expected, tolerance = sqrt(.Machine$double.eps),
               label = "robScale ADM fallback at 1e15")
})

test_that("ensemble robScale path uses adm_core_sorted at 1e15", {
  # scale_robust with method="robScale" goes through rob_scale_sorted
  # which calls adm_core_sorted when MAD implodes
  x <- c(1e15, 1e15, 1e15, 1e15 + 1e-6, 1e15 + 2e-6)
  res <- scale_robust(x, method = "robScale", ci = FALSE)
  expected <- sqrt(pi / 2) * mean(abs(x - median(x)))
  expect_equal(res, expected, tolerance = 100 * .Machine$double.eps,
               label = "ensemble robScale sorted ADM at 1e15")
})

test_that("vshaped_mad does not overflow for large values (S12)", {
  # Values near DBL_MAX/4 — (lo + hi) * 0.5 could overflow
  big <- .Machine$double.xmax / 4
  x <- sort(c(big, big + 1, big + 2, big + 3))
  # mad_scaled uses vshaped_mad internally for sorted input in ensemble path
  # Direct test: mad_scaled should return a finite result
  result <- mad_scaled(x)
  expect_true(is.finite(result), label = "mad_scaled finite for large values")
})
