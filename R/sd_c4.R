#' Unbiased Standard Deviation
#'
#' Computes the sample standard deviation corrected by the \eqn{c_4(n)}{c4(n)}
#' factor to remove the small-sample bias of the square-root estimator.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#'
#' @details
#' The standard \code{\link[stats]{sd}} function computes
#' \eqn{s = \sqrt{s^2}}{s = sqrt(s^2)} where \eqn{s^2} is the unbiased
#' sample variance. However, \eqn{E[s] \neq \sigma}{E[s] != sigma} for finite
#' \eqn{n}; the square root introduces a downward bias.
#'
#' The \eqn{c_4(n)}{c4(n)} correction factor removes this bias:
#'
#' \deqn{\hat\sigma = \frac{s}{c_4(n)} = \frac{s}{\sqrt{2/(n-1)}
#'   \cdot \Gamma(n/2) / \Gamma((n-1)/2)}}{sigma_hat = s / c4(n)}
#'
#' \strong{Statistical Properties.}
#' This is a classical (non-robust) estimator with \bold{100\% ARE} by
#' construction, but \bold{0\% breakdown point} --- a single outlier can make
#' it arbitrarily large.
#'
#' \strong{Numerical Stability.}
#' Uses Welford's online algorithm for numerically stable variance computation,
#' avoiding catastrophic cancellation that affects naive two-pass formulas.
#'
#' @return A single numeric value: the bias-corrected standard deviation.
#'   Returns \code{NA} if \code{n < 2}.
#'
#' @seealso
#' \code{\link[stats]{sd}} for the (biased) sample standard deviation;
#' \code{\link{scale_robust}} for the unified dispatcher;
#' \code{\link{robScale}} for the robust M-estimate of scale.
#'
#' @examples
#' sd_c4(c(1:9))
#'
#' # Compare with base sd() --- difference is small for large n
#' x <- rnorm(1000)
#' c(sd_c4 = sd_c4(x), sd = sd(x))
#'
#' @keywords univar
#' @export
sd_c4 <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (length(x) == 0L) return(NA_real_)
  sd_c4_impl(x)
}
