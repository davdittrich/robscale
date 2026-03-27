# test-robScale-nr-production.R — Production regression pins for NR robScale().
#
# WU-PROD-4: oracle (rob_scale_nr_impl) replaced with pinned numeric values
# recorded after WU-PROD-2 installed NR as the production algorithm.
#
# Pin values (NR era, set.seed as noted):
#   pin1 = 1.1198814880955629825  (seed=42,  n=20)
#   pin2 = 0.57252839884492578992 (x2 worst-case 4-element)
#   pin3 = 0.7079841966702126177  (seed=837, n=4)
#
# Worst-case input (diff ≈ 5.97×10⁻⁸ Aitken vs NR):
#   c(2.048861, -0.954442, 0.28295, 0.20733)
# Empirically found: n=4 inputs with 2 sign changes maximise Aitken/NR divergence.

test_that("Test 1: robScale tight match to NR pin on worst-case 4-element input", {
  x   <- c(2.048861, -0.954442, 0.28295, 0.20733)
  ref <- robScale(x)
  expect_equal(ref, 0.57252839884492578992, tolerance = 1e-9)
})

test_that("Test 2: robScale tight match to NR pin on random seed-837 4-element input", {
  set.seed(837); x <- rnorm(4)
  ref <- robScale(x)
  expect_equal(ref, 0.7079841966702126177, tolerance = 1e-9)
})

test_that("Test 3: robScale matches NR on seed-42 n=20 input (regression guard)", {
  set.seed(42); x <- rnorm(20)
  ref <- robScale(x)
  expect_equal(ref, 1.1198814880955629825, tolerance = 1e-9)
})

test_that("Test 4: robScale fallback (adm_core) still fires for near-degenerate input", {
  # MAD of c(5,5,5,5,6) is 0 → s_init = 0 → adm_core fallback invoked.
  # This path bypasses nr_scale_compute entirely; behaviour must be unchanged.
  result <- robScale(c(5, 5, 5, 5, 6))
  expect_true(length(result) == 1L && is.finite(result) && result > 0)
})

test_that("Test 5: robScale value pin — seed-42 n=20", {
  set.seed(42); x <- rnorm(20)
  ref <- robScale(x)
  expect_equal(ref, 1.1198814880955629825, tolerance = 1e-9)
})
