#' Robust Estimator of Scale Qn
#'
#' Computes the robust estimator of scale \eqn{Q_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x A numeric vector of observations.
#' @param constant Consistency constant. Default is \eqn{2.2191}.
#' @param finite.corr Logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are removed before computation.
#'
#' @return The \eqn{Q_n} estimator of scale.
#'
#' @details
#' \eqn{Q_n} is defined as the \eqn{k}-th order statistic of the \eqn{\binom{n}{2}} absolute
#' differences \eqn{|x_i - x_j|} for \eqn{i < j}, where \eqn{k \approx \binom{n}{2} / 4}.
#' The estimator achieves a 50\% breakdown point and higher statistical efficiency than the
#' Median Absolute Deviation (MAD), without requiring a prior location estimate.
#'
#' The \code{robscale} implementation utilizes a tiered strategy for optimal performance:
#' sorting networks for \eqn{n \le 8}, a specialized Johnson–Mizoguchi algorithm for medium 
#' samples, and cache-aware parallelization for large-scale datasets.
#'
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}(424), 1273--1283.
#' \doi{10.1080/01621459.1993.10476408}
#'
#' @examples
#' qn(c(1:9))
#' x <- c(1, 2, 3, 5, 7, 8)
#' qn(x)
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
