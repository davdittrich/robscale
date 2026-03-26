test_that("scale_robust dispatches to individual estimators exactly", {
  set.seed(42)
  x <- rnorm(10)

  expect_identical(scale_robust(x, method = "gmd"), gmd(x))
  expect_identical(scale_robust(x, method = "sd"),  sd_c4(x))
  expect_identical(scale_robust(x, method = "mad"), mad_scaled(x))
  expect_identical(scale_robust(x, method = "iqr"), iqr_scaled(x))
  expect_identical(scale_robust(x, method = "sn"),  sn(x))
  expect_identical(scale_robust(x, method = "qn"),  qn(x))
  expect_identical(scale_robust(x, method = "robScale"), robScale(x))
})

test_that("scale_robust auto_switch does not intercept named methods (Bug #1)", {
  # With n >= threshold, named methods must return their own value, not gmd.
  set.seed(42)
  x <- rnorm(25)   # n = 25 >= threshold = 20

  for (m in c("sd", "mad", "iqr", "sn", "qn", "robScale")) {
    val_named <- scale_robust(x, method = m, auto_switch = TRUE)
    val_gmd   <- scale_robust(x, method = "gmd", auto_switch = FALSE)
    expect_false(identical(val_named, val_gmd),
                 label = paste0("method=", m, " must not silently return gmd"))
  }

  # gmd with auto_switch = TRUE at n >= threshold should equal gmd(x)
  expect_identical(scale_robust(x, method = "gmd", auto_switch = TRUE), gmd(x))
})

test_that("scale_robust auto-switches ensemble to GMD at threshold", {
  set.seed(42)
  x <- rnorm(25)
  expect_identical(scale_robust(x, threshold = 20L), gmd(x))

  # Below threshold: ensemble is used (different from gmd)
  x_small <- x[1:10]
  expect_false(identical(scale_robust(x_small, threshold = 20L),
                         gmd(x_small)))
})

test_that("scale_robust auto_switch=FALSE forces ensemble for large n", {
  set.seed(42)
  x <- rnorm(50)
  ens     <- scale_robust(x, auto_switch = FALSE)
  gmd_val <- gmd(x)
  expect_false(identical(ens, gmd_val))
})

test_that("scale_robust handles edge cases", {
  expect_true(is.na(scale_robust(numeric(0))))
  expect_true(is.na(scale_robust(1)))
})

test_that("scale_robust handles NA correctly", {
  x <- c(1, 2, 4, NA, 16)
  expect_error(scale_robust(x), "NAs")
  expect_true(is.numeric(scale_robust(x, na.rm = TRUE)))
})

test_that("scale_robust estimates sigma under normality for all methods", {
  set.seed(42)
  y <- rnorm(5000)
  for (m in c("gmd", "sd", "mad", "iqr", "sn", "qn", "robScale")) {
    val <- scale_robust(y, method = m, auto_switch = FALSE)
    expect_equal(val, 1.0, tolerance = 0.1, label = m)
  }
})

test_that("scale_robust ci = TRUE analytical path returns robscale_ci", {
  set.seed(42)
  x <- rnorm(15)
  for (m in c("gmd", "sd", "mad", "iqr", "sn", "qn", "robScale")) {
    r <- scale_robust(x, method = m, ci = TRUE, boot_method = "analytical")
    expect_s3_class(r, "robscale_ci")
    expect_equal(r$boot_method, "analytical")
    expect_true(is.finite(r$ci["lower"]))
    expect_true(is.finite(r$ci["upper"]))
    expect_true(r$ci["lower"] <= r$ci["upper"])
  }
})

test_that("scale_robust ci = TRUE bootstrap path returns robscale_ci", {
  set.seed(42)
  x <- rnorm(15)
  for (m in c("gmd", "sd", "mad", "iqr", "sn", "qn", "robScale")) {
    r <- scale_robust(x, method = m, ci = TRUE, boot_method = "percentile",
                      n_boot = 100L)
    expect_s3_class(r, "robscale_ci")
    expect_equal(r$boot_method, "percentile")
    expect_true(is.finite(r$ci["lower"]))
    expect_true(is.finite(r$ci["upper"]))
    expect_true(r$ci["lower"] <= r$ci["upper"])
  }
})

test_that("scale_robust method name normalisation in CI output", {
  set.seed(42)
  x <- rnorm(15)
  expect_equal(scale_robust(x, method = "sd",  ci = TRUE)$method, "sd_c4")
  expect_equal(scale_robust(x, method = "mad", ci = TRUE)$method, "mad_scaled")
  expect_equal(scale_robust(x, method = "iqr", ci = TRUE)$method, "iqr_scaled")
  expect_equal(scale_robust(x, method = "gmd", ci = TRUE)$method, "gmd")
  expect_equal(scale_robust(x, method = "sn",  ci = TRUE)$method, "sn")
})

test_that("scale_robust: analytical CI not available for ensemble", {
  set.seed(42)
  x <- rnorm(10)
  expect_error(
    scale_robust(x, method = "ensemble", ci = TRUE,
                 boot_method = "analytical"),
    "analytical.*ensemble"
  )
})

test_that("scale_robust ci = TRUE auto-switched to gmd still returns robscale_ci", {
  set.seed(42)
  x <- rnorm(25)
  r <- scale_robust(x, method = "ensemble", auto_switch = TRUE, ci = TRUE)
  expect_s3_class(r, "robscale_ci")
  expect_equal(r$method, "gmd")
})

test_that("get_consistency_constant returns correct values", {
  expect_equal(get_consistency_constant("mad"), 1.482602218505602)
  expect_equal(get_consistency_constant("gmd"), 0.886226925452758)
  expect_equal(get_consistency_constant("iqr"), 0.741301109252801)
})

test_that("get_consistency_constant requires n for sample-dependent methods", {
  expect_error(get_consistency_constant("c4"))
  expect_error(get_consistency_constant("sn"))
  expect_error(get_consistency_constant("qn"))
  expect_true(is.numeric(get_consistency_constant("c4", 5)))
  expect_true(is.numeric(get_consistency_constant("sn", 10)))
  expect_true(is.numeric(get_consistency_constant("qn", 10)))
})
