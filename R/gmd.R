#' Gini Mean Difference
#'
#' Computes the Gini mean difference, scaled by a consistency constant for
#' asymptotic normality under the Gaussian model.
#'
#' @param x A numeric vector.
#' @param constant Consistency constant for asymptotic normality at the
#'   Gaussian. Defaults to \eqn{\sqrt{\pi}/2 \approx 0.8862}{sqrt(pi)/2 ~
#'   0.8862}. Set to \code{1} for the raw (unscaled) Gini mean difference.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @details
#' The Gini mean difference is defined as
#'
#' \deqn{\mathrm{GMD}(x) = C \cdot \frac{2}{n(n-1)}
#'   \sum_{i=1}^{n} (2i - n - 1)\, x_{(i)}}{GMD(x) = C * 2/(n(n-1)) * sum((2i - n - 1) * x_(i))}
#'
#' where \eqn{x_{(1)} \le \ldots \le x_{(n)}} are the order statistics and
#' \eqn{C} is the consistency constant. The computation requires a full sort
#' (O(n log n)).
#'
#' \strong{Statistical Properties.}
#' The default constant \eqn{C = \sqrt{\pi}/2}{C = sqrt(pi)/2} ensures that the
#' GMD is a consistent estimator of \eqn{\sigma} under the Gaussian model. The
#' GMD achieves an \bold{asymptotic relative efficiency (ARE) of 0.98} compared
#' to the sample standard deviation, making it the most efficient robust
#' alternative in this package. Its breakdown point is \eqn{1 - 1/\sqrt{2}}{1 - 1/sqrt(2)}, approximately 29.3\%.
#'
#' @return If \code{ci = FALSE} (default), a single numeric value: the scaled
#'   Gini mean difference. Returns \code{0} if \code{n < 2}; returns \code{NA}
#'   if \code{x} has length zero after removal of \code{NA}s.
#'   If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
#'
#' @references
#' Gini, C. (1912) \emph{Variabilita e mutabilita}. Bologna.
#'
#' @seealso
#' \code{\link{scale_robust}} for the unified dispatcher;
#' \code{\link{qn}} and \code{\link{sn}} for high-breakdown scale estimators.
#'
#' @examples
#' gmd(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' gmd(x)                      # with consistency constant
#' gmd(x, constant = 1)        # raw Gini mean difference
#'
#' # Asymptotic confidence interval
#' gmd(x, ci = TRUE)
#'
#' @keywords univar robust
#' @export
gmd <- function(x, constant = 0.886226925452758, na.rm = FALSE,
                ci = FALSE, level = 0.95) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }
  res <- .Call(`_robscale_gmd_impl`, x, constant)
  if (ci) return(.analytical_ci(res, n, are = .are_values[["gmd"]], level, "gmd"))
  res
}
