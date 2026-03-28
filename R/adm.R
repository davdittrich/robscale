#' Average Distance to the Median
#'
#' Computes the mean absolute deviation from the median, scaled by a consistency
#' constant for asymptotic normality under the Gaussian model.
#'
#' @param x A numeric vector.
#' @param center Optional numeric scalar giving the central value from which to
#'   measure the average absolute distance.  Defaults to the median of
#'   \code{x}.
#' @param constant Consistency constant for asymptotic normality at the
#'   Gaussian.  Defaults to \eqn{\sqrt{\pi/2} \approx 1.2533}{sqrt(pi/2) ~
#'   1.2533} (Nair, 1947).  Set to \code{1} for the raw (unscaled) mean
#'   absolute deviation.
#' @param na.rm Logical.  If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation.  If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @details
#' The average distance to the median (ADM) is defined as
#'
#' \deqn{\mathrm{ADM}(x) \;=\; C \cdot \frac{1}{n}\sum_{i=1}^{n}
#'   |x_i - \mathrm{med}(x)|}{ADM(x) = C * mean(|x - med(x)|)}
#'
#' where \eqn{C} is the consistency constant and \eqn{\mathrm{med}(x)}{median(x)}
#' is the sample median. When \code{center} is supplied, it replaces the
#' sample median.
#'
#' \strong{Statistical Properties.}
#' The default constant \eqn{C = \sqrt{\pi/2}}{C = sqrt(pi/2)} ensures that the
#' ADM is a consistent estimator of the standard deviation \eqn{\sigma} under
#' the Gaussian model. At the normal distribution, the ADM achieves an
#' \bold{asymptotic relative efficiency (ARE) of 0.88} compared to the sample
#' standard deviation.
#'
#' While the ADM is less efficient than the standard deviation for purely
#' Gaussian data, it offers superior resistance to "implosion" (the estimate
#' collapsing to zero). Its implosion breakdown point is
#' \eqn{(n-1)/n}{(n - 1) / n}, meaning it only collapses if all but one
#' observation are identical. Conversely, its explosion breakdown point is
#' \eqn{1/n}, similar to the sample mean. These properties make the ADM the
#' ideal \bold{implosion fallback} for the M-estimator of scale in
#' \code{\link{robScale}}.
#'
#' \strong{Computational Performance.}
#' This implementation employs a tiered selection strategy: optimal sorting
#' networks for \eqn{n \le 16} and adaptive \eqn{O(n)} selection
#' (Floyd--Rivest or pdqselect, depending on cache-derived thresholds) for
#' larger datasets. This avoids the full sort required by the standard
#' \code{\link{median}} function.
#'
#' @return A single numeric value: the scaled mean absolute deviation from the
#'   center.  Returns \code{NA} if \code{x} has length zero after removal of
#'   \code{NA}s.
#'
#' @references
#' Nair, K. R. (1947) A Note on the Mean Deviation from the Median.
#' \emph{Biometrika}, \bold{34}(3/4), 360--362. \doi{10.2307/2332448}
#'
#' Rousseeuw, P. J. and Verboven, S. (2002) Robust estimation in very small
#' samples. \emph{Computational Statistics & Data Analysis}, \bold{40}(4),
#' 741--758. \doi{10.1016/S0167-9473(02)00078-6}
#'
#' @seealso
#' \code{\link[stats]{mad}} for the median absolute deviation from the
#' \code{\link{median}};
#' \code{\link{robScale}} for the M-estimator of scale that uses the ADM as
#' an implosion fallback.
#'
#' @examples
#' adm(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' adm(x)                      # with consistency constant
#' adm(x, constant = 1)        # raw mean absolute deviation
#'
#' # Supply a known center
#' adm(x, center = 4.0)
#'
#' @keywords univar robust
#' @export
adm <- function(x, center = NULL, constant = 1.2533141373155001, na.rm = FALSE,
                ci = FALSE, level = 0.95) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (any(!is.finite(x))) stop("'x' must not contain non-finite values (Inf, -Inf, NaN)")
  if (length(x) == 0L) return(NA_real_)
  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }
  if (is.null(center)) {
    res <- .Call(`_robscale_adm_impl_auto`, x, constant)
  } else {
    res <- .Call(`_robscale_adm_impl`, x, center, constant)
  }
  if (ci) return(.analytical_ci(res, length(x), are = .are_values[["adm"]], level, "adm"))
  res
}
