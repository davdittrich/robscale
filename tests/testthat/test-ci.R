# --- Individual estimator CI tests ---

test_that("each estimator with ci=TRUE returns robscale_ci class", {
  set.seed(42)
  x <- rnorm(50)
  fns <- list(sn = sn, qn = qn, gmd = gmd, mad_scaled = mad_scaled,
              iqr_scaled = iqr_scaled, sd_c4 = sd_c4, robScale = robScale)
  for (nm in names(fns)) {
    res <- fns[[nm]](x, ci = TRUE)
    expect_s3_class(res, "robscale_ci")
    expect_true(is.numeric(res$estimate), label = nm)
    expect_named(res$ci, c("lower", "upper"), label = nm)
    expect_equal(res$level, 0.95, label = nm)
  }
})

test_that("CI contains the point estimate", {
  set.seed(42)
  x <- rnorm(50)
  fns <- list(sn = sn, qn = qn, gmd = gmd, mad_scaled = mad_scaled,
              iqr_scaled = iqr_scaled, sd_c4 = sd_c4, robScale = robScale)
  for (nm in names(fns)) {
    res <- fns[[nm]](x, ci = TRUE)
    expect_true(res$ci[["lower"]] < res$estimate, label = nm)
    expect_true(res$ci[["upper"]] > res$estimate, label = nm)
  }
})

test_that("level=0.99 gives wider CI than level=0.90", {
  set.seed(42)
  x <- rnorm(50)
  fns <- list(sn = sn, qn = qn, gmd = gmd, mad_scaled = mad_scaled,
              iqr_scaled = iqr_scaled, sd_c4 = sd_c4, robScale = robScale)
  for (nm in names(fns)) {
    ci90 <- fns[[nm]](x, ci = TRUE, level = 0.90)
    ci99 <- fns[[nm]](x, ci = TRUE, level = 0.99)
    w90 <- ci90$ci[["upper"]] - ci90$ci[["lower"]]
    w99 <- ci99$ci[["upper"]] - ci99$ci[["lower"]]
    expect_true(w99 > w90, label = nm)
  }
})

test_that("sd_c4 CI is asymmetric (chi-squared)", {
  set.seed(42)
  x <- rnorm(20)
  res <- sd_c4(x, ci = TRUE)
  lower_dist <- res$estimate - res$ci[["lower"]]
  upper_dist <- res$ci[["upper"]] - res$estimate
  expect_false(isTRUE(all.equal(lower_dist, upper_dist)))
  expect_equal(res$method, "sd_c4")
})

test_that("scalar return unchanged when ci=FALSE", {
  set.seed(42)
  x <- rnorm(50)
  expect_true(is.numeric(sn(x)))
  expect_true(is.atomic(sn(x)))
  expect_length(sn(x), 1)
  expect_true(is.numeric(sd_c4(x)))
  expect_true(is.atomic(sd_c4(x)))
})

test_that("edge case: n=2 produces valid CIs", {
  x <- c(1.0, 3.0)
  fns <- list(sn = sn, qn = qn, gmd = gmd, mad_scaled = mad_scaled,
              iqr_scaled = iqr_scaled, sd_c4 = sd_c4)
  for (nm in names(fns)) {
    res <- fns[[nm]](x, ci = TRUE)
    expect_s3_class(res, "robscale_ci")
    expect_true(is.finite(res$ci[["lower"]]), label = nm)
    expect_true(is.finite(res$ci[["upper"]]), label = nm)
    expect_true(res$ci[["lower"]] < res$ci[["upper"]], label = nm)
  }
})

test_that("print.robscale_ci produces output", {
  set.seed(42)
  x <- rnorm(50)
  res <- sn(x, ci = TRUE)
  expect_output(print(res), "sn estimate:")
  expect_output(print(res), "95% CI")
})

test_that("analytical CI method names are correct", {
  set.seed(42)
  x <- rnorm(50)
  expect_equal(sn(x, ci = TRUE)$method, "sn")
  expect_equal(qn(x, ci = TRUE)$method, "qn")
  expect_equal(gmd(x, ci = TRUE)$method, "gmd")
  expect_equal(mad_scaled(x, ci = TRUE)$method, "mad_scaled")
  expect_equal(iqr_scaled(x, ci = TRUE)$method, "iqr_scaled")
  expect_equal(sd_c4(x, ci = TRUE)$method, "sd_c4")
  expect_equal(robScale(x, ci = TRUE)$method, "robScale")
})
