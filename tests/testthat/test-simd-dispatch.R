# --- Config assertions ---

test_that("has_avx2 is present in get_qnsn_config()", {
  config <- robscale:::get_qnsn_config()
  expect_true("has_avx2" %in% names(config))
  expect_type(config$has_avx2, "logical")
})

test_that("has_avx2 is consistent with simd_level", {
  config <- robscale:::get_qnsn_config()
  arch <- Sys.info()[["machine"]]
  if (arch %in% c("arm64", "aarch64")) {
    expect_false(config$has_avx2)
  } else if (arch %in% c("x86_64", "amd64")) {
    expect_equal(config$has_avx2, config$simd_level >= 2L)
  }
})

# --- Correctness under runtime dispatch ---

test_that("robScale is valid across sizes (exercises bulk_tanh dispatch)", {
  set.seed(123)
  for (n in c(100, 1000, 10000)) {
    x <- rnorm(n)
    result <- robScale(x)
    expect_true(is.finite(result), label = paste("robScale finite at n =", n))
    expect_true(result > 0, label = paste("robScale positive at n =", n))
  }
})

test_that("qn is valid across brute-force sizes (exercises AVX2 kernel dispatch)", {
  set.seed(123)
  for (n in c(17, 32, 48, 64)) {
    x <- rnorm(n)
    result <- qn(x)
    expect_true(is.finite(result), label = paste("qn finite at n =", n))
    expect_true(result > 0, label = paste("qn positive at n =", n))
  }
})

test_that("robLoc is valid across sizes (exercises bulk_tanh dispatch)", {
  set.seed(123)
  for (n in c(100, 1000)) {
    x <- rnorm(n)
    result <- robLoc(x)
    expect_true(is.finite(result), label = paste("robLoc finite at n =", n))
  }
})

# --- Determinism (no path-dependent numerical drift) ---

test_that("robScale is deterministic across repeated calls", {
  set.seed(42)
  x <- rnorm(5000)
  results <- replicate(50, robScale(x))
  expect_equal(length(unique(results)), 1L)
})

test_that("qn is deterministic across repeated calls", {
  set.seed(42)
  x <- rnorm(64)
  results <- replicate(50, qn(x))
  expect_equal(length(unique(results)), 1L)
})
