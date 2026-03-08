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
  n <- length(x)
  if (n < 2) return(NA_real_)
  
  if (is.double(x)) {
    if (constant == 1.1926 && finite.corr) return(.Call(`_robscale_C_sn_fast`, x))
    res <- .Call(`_robscale_C_sn_fast`, x)
  } else if (is.integer(x)) {
    if (constant == 1.1926 && finite.corr) return(.Call(`_robscale_C_sn_int_fast`, x))
    res <- .Call(`_robscale_C_sn_int_fast`, x)
  } else {
    x <- as.double(x)
    if (constant == 1.1926 && finite.corr) return(.Call(`_robscale_C_sn_fast`, x))
    res <- .Call(`_robscale_C_sn_fast`, x)
  }
  
  if (constant != 1.1926) {
    res <- res * (constant / 1.19259855312321)
  }
  
  if (!finite.corr) {
    res <- res / .Call(`_robscale_C_get_sn_factor`, n)
  }
  
  res
}
