test_that("scale_robust dispatches to individual estimators exactly", {
  set.seed(42)
  x <- rnorm(10)

  expect_identical(scale_robust(x, method = "gmd", auto_switch = FALSE), gmd(x))
  expect_identical(scale_robust(x, method = "sd", auto_switch = FALSE), sd_c4(x))
  expect_identical(scale_robust(x, method = "mad", auto_switch = FALSE), mad_scaled(x))
  expect_identical(scale_robust(x, method = "iqr", auto_switch = FALSE), iqr_scaled(x))
  expect_identical(scale_robust(x, method = "sn", auto_switch = FALSE), sn(x))
  expect_identical(scale_robust(x, method = "qn", auto_switch = FALSE), qn(x))
  expect_identical(scale_robust(x, method = "robScale", auto_switch = FALSE), robScale(x))
})

test_that("scale_robust auto-switches to GMD at threshold", {
  set.seed(42)
  x <- rnorm(25)
  expect_identical(scale_robust(x, threshold = 20L), gmd(x))
  # Below threshold, should use ensemble (different from gmd)
  x_small <- x[1:10]
  expect_false(identical(scale_robust(x_small, threshold = 20L),
                         gmd(x_small)))
})

test_that("scale_robust auto_switch=FALSE forces method", {
  set.seed(42)
  x <- rnorm(50)
  # Without auto_switch, ensemble is used even for large n
  ens <- scale_robust(x, auto_switch = FALSE)
  gmd_val <- gmd(x)
  # Ensemble != GMD (they differ for large n)
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
  methods <- c("gmd", "sd", "mad", "iqr", "sn", "qn", "robScale")
  for (m in methods) {
    val <- scale_robust(y, method = m, auto_switch = FALSE)
    expect_equal(val, 1.0, tolerance = 0.1, label = m)
  }
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
