# WU-7: R input validation + diag.cpp NR update
# Covers: R6 (sd_c4 n=1), R7 (qn/sn NA), R8 (is.numeric), R9 (diag NR),
#         S6 (fallback), S7 (level), S8 (Inf), S9 (scalar params)

# --- R6: sd_c4 n=1 guard ---
test_that("sd_c4 returns NA_real_ for n=1 (R6)", {
  expect_identical(sd_c4(42), NA_real_)
})

test_that("sd_c4 returns NA_real_ for n=0 (R6 regression)", {
  expect_identical(sd_c4(numeric(0)), NA_real_)
})

# --- R7: qn/sn NA handling ---
test_that("qn errors on NA with na.rm=FALSE (R7)", {
  expect_error(qn(c(1, 2, NA)), "NAs")
})

test_that("sn errors on NA with na.rm=FALSE (R7)", {
  expect_error(sn(c(1, 2, NA)), "NAs")
})

test_that("qn works with na.rm=TRUE (R7 regression)", {
  expect_true(is.finite(qn(c(1, 2, NA, 3, 4), na.rm = TRUE)))
})

test_that("sn works with na.rm=TRUE (R7 regression)", {
  expect_true(is.finite(sn(c(1, 2, NA, 3, 4), na.rm = TRUE)))
})

# --- R8: is.numeric validation ---
test_that("all wrappers error on logical input (R8)", {
  fns <- list(sd_c4 = sd_c4, qn = qn, sn = sn, gmd = gmd,
              mad_scaled = mad_scaled, adm = adm, robScale = robScale,
              robLoc = robLoc, iqr_scaled = iqr_scaled)
  for (nm in names(fns)) {
    expect_error(fns[[nm]](c(TRUE, FALSE, TRUE, TRUE, FALSE)),
                 "numeric", label = paste0(nm, " rejects logical"))
  }
})

test_that("robScale errors on character input (R8)", {
  expect_error(robScale(letters[1:5]), "numeric")
})

test_that("robScale accepts integer input (R8 regression)", {
  expect_true(is.finite(robScale(1:10)))
})

# --- R9: diag.cpp NR (not Aitken) ---
test_that("rob_scale_diag_impl uses NR and matches robScale (R9)", {
  for (n in c(4, 10, 20, 50, 100)) {
    set.seed(42 + n)
    x <- rnorm(n)
    diag_res <- rob_scale_diag_impl(x)
    prod_res <- robScale(x)
    expect_equal(diag_res$scale, prod_res,
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste0("diag matches robScale at n=", n))
    expect_equal(diag_res$aitken_fires, 0L,
                 label = paste0("no Aitken at n=", n))
    expect_lte(diag_res$outer_iters, 15L,
               label = paste0("NR converges in <=15 iters at n=", n))
  }
})

# --- S6: fallback simplification ---
test_that("robScale fallback='adm' works (S6)", {
  x <- c(5, 5, 5, 5, 5)
  expect_equal(robScale(x, fallback = "adm"), 0)
})

test_that("robScale fallback='na' works (S6)", {
  x <- c(5, 5, 5, 5, 5)
  expect_true(is.na(robScale(x, fallback = "na")))
})

test_that("robScale default fallback works (S6 regression)", {
  expect_true(is.finite(robScale(rnorm(10))))
})

# --- S7: level validation ---
test_that("CI wrappers error on invalid level (S7)", {
  x <- 1:10
  expect_error(sd_c4(x, ci = TRUE, level = 0), "level")
  expect_error(sd_c4(x, ci = TRUE, level = 1.1), "level")
  expect_error(robScale(x, ci = TRUE, level = -0.5), "level")
  expect_error(qn(x, ci = TRUE, level = 2), "level")
})

# --- S8: Inf/-Inf guard ---
test_that("wrappers error on Inf input (S8)", {
  expect_error(robScale(c(1, 2, Inf, 3, 4)), "finite")
  expect_error(sd_c4(c(1, -Inf, 3, 4, 5)), "finite")
  expect_error(gmd(c(1, 2, Inf)), "finite")
})

# --- S9: vector-valued scalar params ---
test_that("robScale errors on vector loc (S9)", {
  expect_error(robScale(1:5, loc = c(2, 3)), "single")
})
