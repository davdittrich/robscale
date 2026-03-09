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
#' The \eqn{Q_n} estimator is defined as the \eqn{k}-th order statistic of the
#' \eqn{\binom{n}{2}} absolute pairwise differences \eqn{|x_i - x_j|} for
#' \eqn{i < j}. Specifically, \eqn{Q_n = C \cdot d_{(k)}} where \eqn{d} is the
#' set of absolute differences and \eqn{k = \binom{h}{2}} with
#' \eqn{h = \lfloor n/2 \rfloor + 1}.
#'
#' \strong{Statistical Properties.}
#' \eqn{Q_n} is a highly robust estimator with a \bold{50\% breakdown point}.
#' Unlike the Median Absolute Deviation (MAD), \eqn{Q_n} does not require a
#' prior location estimate, making it more robust in asymmetric distributions.
#' At the Gaussian distribution, it achieves an \bold{asymptotic relative
#' efficiency (ARE) of 0.82}, significantly higher than the 0.37 achieved by the
#' MAD.
#'
#' \strong{Computational Performance.}
#' While the naive calculation of \eqn{Q_n} requires \eqn{O(n^2)} space and
#' time, this implementation employs a specialized \eqn{O(n \log n)} algorithm.
#' The package utilizes a tiered execution strategy:
#' \itemize{
#'   \item \bold{Optimal sorting networks} for very small samples (\eqn{n \le 8}).
#'     These networks eliminate branch misprediction in the target regimes of
#'     extremely small samples.
#'   \item A specialized \bold{Johnson--Mizoguchi selection} algorithm for
#'     medium-sized datasets.
#'   \item \bold{Cache-aware parallelization} via Intel TBB (Threading Building
#'     Blocks, when \code{ROBSCALE_FAST=1}) for large-scale data.
#' }
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
