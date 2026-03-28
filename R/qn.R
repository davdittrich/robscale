#' Robust Estimator of Scale Qn
#'
#' Computes the robust estimator of scale \eqn{Q_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x A numeric vector of observations.
#' @param constant Consistency constant. Default is \eqn{2.2191}.
#' @param finite.corr Logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are removed before computation.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @return If \code{ci = FALSE} (default), the \eqn{Q_n} estimator of scale
#'   (a scalar). If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
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
#' The implementation uses a tiered execution strategy:
#' \itemize{
#'   \item \bold{Optimal sorting networks} for very small samples (\eqn{n \le 16}).
#'     These networks eliminate branch misprediction in the target regimes of
#'     extremely small samples.
#'   \item A \bold{Croux--Rousseeuw weighted-median refinement} algorithm
#'     for medium and large datasets.
#'   \item \bold{Cache-aware parallelization} via Intel TBB (Threading Building
#'     Blocks) for large-scale data, with thresholds derived from the detected
#'     per-core L2 cache size.
#' }
#'
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}(424), 1273--1283.
#' \doi{10.1080/01621459.1993.10476408}
#'
#' Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators.
#' \emph{arXiv preprint arXiv:2209.12268}.
#'
#' @seealso \code{\link{sn}} for the \eqn{S_n} scale estimator;
#'   \code{\link{robScale}} for the M-estimator of scale;
#'   \code{\link{adm}} for the average distance to median.
#'
#' @examples
#' qn(c(1:9))
#' x <- c(1, 2, 3, 5, 7, 8)
#' qn(x)
#'
#' # Asymptotic confidence interval
#' qn(x, ci = TRUE)
#'
#' @export
qn <- function(x, constant = 2.21914446598508, finite.corr = TRUE, na.rm = FALSE,
               ci = FALSE, level = 0.95) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)

  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }

  fast <- !ci && constant == 2.21914446598508 && finite.corr

  res <- if (is.double(x)) .Call(`_robscale_C_qn_fast`, x)
         else .Call(`_robscale_C_qn_int_fast`, x)
  if (fast) return(res)

  if (constant != 2.21914446598508) {
    res <- res * (constant / 2.21914446598508)
  }

  if (!finite.corr) {
    res <- res / .Call(`_robscale_C_get_qn_factor`, n)
  }

  if (ci) return(.analytical_ci(res, n, are = .are_values[["qn"]], level, "qn"))
  res
}
