#' Robust Estimator of Scale Qn
#'
#' Computes the robust estimator of scale \eqn{Q_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x a numeric vector of observations.
#' @param constant consistency constant. Default is \eqn{2.2191}.
#' @param finite.corr logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm logical; if \code{TRUE}, \code{NA} values are removed before computation.
#'
#' @return The \eqn{Q_n} estimator of scale.
#'
#' @details
#' \eqn{Q_n} is defined as \eqn{2.2191 \cdot \{ |x_i - x_j|; i < j \}_{(k)}} where \eqn{k \approx \binom{n}{2} / 4}.
#' It has a breakdown point of 50%, is more efficient than the MAD, and does not require a location estimate.
#'
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#'
#' @export
qn <- function(x, constant = 2.2191, finite.corr = TRUE, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)
  
  if (is.double(x)) {
    if (identical(constant, 2.2191) && finite.corr) return(.Call(`_robscale_C_qn_fast`, x))
    res <- .Call(`_robscale_C_qn_fast`, x)
  } else if (is.integer(x)) {
    if (identical(constant, 2.2191) && finite.corr) return(.Call(`_robscale_C_qn_int_fast`, x))
    res <- .Call(`_robscale_C_qn_int_fast`, x)
  } else {
    x <- as.double(x)
    if (identical(constant, 2.2191) && finite.corr) return(.Call(`_robscale_C_qn_fast`, x))
    res <- .Call(`_robscale_C_qn_fast`, x)
  }
  
  if (!identical(constant, 2.2191)) {
    res <- res * (constant / 2.21914446598508)
  }
  
  if (!finite.corr) {
    res <- res / .Call(`_robscale_C_get_qn_factor`, n)
  }
  
  res
}
