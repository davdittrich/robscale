#' Robust Estimator of Scale Sn
#'
#' Computes the robust estimator of scale \eqn{S_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x A numeric vector of observations.
#' @param constant Consistency constant. Default is \eqn{1.1926}.
#' @param finite.corr Logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are removed before computation.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @return If \code{ci = FALSE} (default), the \eqn{S_n} estimator of scale
#'   (a scalar). If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
#'
#' @details
#' The \eqn{S_n} estimator is defined as the median of medians:
#'
#' \deqn{S_n = C \cdot \text{med}_i \left\{ \text{med}_j |x_i - x_j| \right\}}{Sn = C * med_i(med_j(|x_i - x_j|))}
#'
#' \strong{Statistical Properties.}
#' \eqn{S_n} achieves a \bold{50\% breakdown point} and provides superior
#' statistical efficiency compared to the Median Absolute Deviation (MAD). At
#' the Gaussian distribution, it has an \bold{asymptotic relative efficiency
#' (ARE) of 0.58}, which is significantly higher than the 0.37 of the MAD.
#' Unlike M-estimators of scale, \eqn{S_n} is an explicit function of the
#' data and does not require an iterative solution or a prior location
#' estimate.
#'
#' \strong{Computational Performance.}
#' This implementation provides a highly optimized \bold{\eqn{O(n \log n)}}
#' algorithm, avoiding the \eqn{O(n^2)} complexity of the naive definition.
#' The execution strategy is tiered for maximum efficiency:
#' \itemize{
#'   \item \bold{Optimal sorting networks} are used for the target regime of
#'     very small samples (\eqn{n \le 16}). These networks minimize control
#'     flow overhead and maximize pipeline utilization.
#'   \item \bold{Highly tuned kernels} are employed for general datasets,
#'     leveraging C++17 features for cache-aware computation.
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
#' @seealso \code{\link{qn}} for the \eqn{Q_n} scale estimator;
#'   \code{\link{robScale}} for the M-estimator of scale;
#'   \code{\link{adm}} for the average distance to median.
#'
#' @examples
#' sn(c(1:9))
#' x <- c(1, 2, 3, 5, 7, 8)
#' sn(x)
#'
#' # Asymptotic confidence interval
#' sn(x, ci = TRUE)
#'
#' @export
sn <- function(x, constant = 1.19259855312321, finite.corr = TRUE, na.rm = FALSE,
               ci = FALSE, level = 0.95) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(NA_real_)

  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }

  fast <- !ci && constant == 1.19259855312321 && finite.corr

  if (is.double(x)) {
    if (fast) return(.Call(`_robscale_C_sn_fast`, x))
    res <- .Call(`_robscale_C_sn_fast`, x)
  } else {
    if (fast) return(.Call(`_robscale_C_sn_int_fast`, x))
    res <- .Call(`_robscale_C_sn_int_fast`, x)
  }

  if (constant != 1.19259855312321) {
    res <- res * (constant / 1.19259855312321)
  }

  if (!finite.corr) {
    res <- res / .Call(`_robscale_C_get_sn_factor`, n)
  }

  if (ci) return(.analytical_ci(res, n, are = .are_values[["sn"]], level, "sn"))
  res
}
