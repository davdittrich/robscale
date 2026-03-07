#' Robust Estimator of Scale Sn
#'
#' Computes the robust estimator of scale \eqn{S_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x a numeric vector of observations.
#' @param constant consistency constant. Default is \eqn{1.1926}.
#' @param finite.corr logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm logical; if \code{TRUE}, \code{NA} values are removed before computation.
#'
#' @return The \eqn{S_n} estimator of scale.
#'
#' @details
#' \eqn{S_n} is defined as \eqn{1.1926 \cdot \text{med}_i \{ \text{med}_j |x_i - x_j| \}}.
#' It has a breakdown point of 50% and is more efficient than the MAD.
#'
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#'
#' @export
sn <- function(x, constant = 1.1926, finite.corr = TRUE, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  
  if (is.double(x)) {
    res <- C_sn_fast(x)
  } else if (is.integer(x)) {
    res <- C_sn_int_fast(x)
  } else {
    res <- C_sn_fast(as.double(x))
  }
  
  # The C++ implementation already applies CONST_SN and get_sn_factor.
  # If the user provided a custom constant, we adjust.
  if (constant != 1.1926) {
    res <- res * (constant / 1.19259855312321)
  }
  
  if (!finite.corr) {
    # If no finite correction, we divide by the factor applied in C++
    n <- length(x)
    res <- res / robscale:::get_sn_factor_r(n)
  }
  
  res
}

# Helper for finite correction (internal)
get_sn_factor_r <- function(n) {
  if (n <= 100) {
    factors <- c(0.00000, 0.74303, 1.84983, 0.95505, 1.34857, 0.99413, 1.19832, 1.00496, 1.13178, 1.00689,
                1.09592, 1.00635, 1.07423, 1.00513, 1.06006, 1.00384, 1.05006, 1.00281, 1.04297, 1.00219,
                1.03738, 1.00139, 1.03311, 1.00091, 1.02969, 1.00066, 1.02686, 1.00045, 1.02449, 1.00005,
                1.02260, 0.99995, 1.02087, 0.99974, 1.01950, 0.99978, 1.01830, 0.99960, 1.01717, 0.99969,
                1.01619, 0.99960, 1.01538, 0.99955, 1.01460, 0.99960, 1.01391, 0.99948, 1.01324, 0.99953,
                1.01264, 0.99954, 1.01228, 0.99949, 1.01175, 0.99950, 1.01127, 0.99955, 1.01090, 0.99959,
                1.01054, 0.99954, 1.01023, 0.99963, 1.00988, 0.99968, 1.00951, 0.99959, 1.00923, 0.99966,
                1.00902, 0.99965, 1.00877, 0.99964, 1.00851, 0.99966, 1.00835, 0.99968, 1.00810, 0.99966,
                1.00790, 0.99970, 1.00765, 0.99970, 1.00762, 0.99968, 1.00740, 0.99972, 1.00723, 0.99973,
                1.00705, 0.99974, 1.00689, 0.99974, 1.00674, 0.99978, 1.00661, 0.99973, 1.00650, 0.99982)
    return(factors[n])
  }
  inv_n <- 1.0 / n
  if (n %% 2 == 1) {
    return(1.0 + 0.7096 * inv_n - 7.3604 * inv_n^2)
  } else {
    return(1.0 + 0.0391 * inv_n - 6.1719 * inv_n^2)
  }
}
