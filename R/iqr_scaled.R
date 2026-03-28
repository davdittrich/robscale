#' Scaled Interquartile Range
#'
#' Computes the interquartile range, scaled by a consistency constant for
#' asymptotic normality under the Gaussian model.
#'
#' @param x A numeric vector.
#' @param constant Consistency constant for asymptotic normality at the
#'   Gaussian. Defaults to
#'   \eqn{1/(\Phi^{-1}(0.75) - \Phi^{-1}(0.25)) \approx 0.7413}{1/(Phi^-1(0.75) - Phi^-1(0.25)) ~ 0.7413}.
#'   Set to \code{1} for the raw interquartile range.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @details
#' The scaled IQR is defined as
#'
#' \deqn{\mathrm{IQR}_s(x) = C \cdot (Q_{0.75} - Q_{0.25})}{IQR_s(x) = C * (Q_0.75 - Q_0.25)}
#'
#' where \eqn{Q_p} denotes the Type 7 quantile (R default) and \eqn{C} is
#' the consistency constant.
#'
#' \strong{Statistical Properties.}
#' The default constant ensures consistency for \eqn{\sigma} under the Gaussian
#' model. The IQR achieves an \bold{asymptotic relative efficiency (ARE) of
#' 0.37} compared to the sample standard deviation. Its breakdown point is 25\%.
#'
#' \strong{Computational Performance.}
#' Unlike \code{\link[stats]{IQR}}, which requires a full sort, this
#' implementation uses O(n) selection via the pdqselect algorithm for each
#' quantile, providing a substantial speedup for large datasets.
#'
#' @return If \code{ci = FALSE} (default), a single numeric value: the scaled
#'   interquartile range. Returns \code{0} if \code{n < 2}; returns \code{NA}
#'   if \code{x} has length zero after removal of \code{NA}s.
#'   If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
#'
#' @seealso
#' \code{\link[stats]{IQR}} for the base R (unscaled) interquartile range;
#' \code{\link{scale_robust}} for the unified dispatcher;
#' \code{\link{mad_scaled}} for the scaled median absolute deviation.
#'
#' @examples
#' iqr_scaled(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' iqr_scaled(x)                 # with consistency constant
#' iqr_scaled(x, constant = 1)   # raw IQR
#'
#' # Asymptotic confidence interval
#' iqr_scaled(x, ci = TRUE)
#'
#' @keywords univar robust
#' @export
iqr_scaled <- function(x, constant = 0.741301109252801, na.rm = FALSE,
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
  n <- length(x)
  if (n == 0L) return(NA_real_)
  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }
  res <- iqr_impl(x, constant)
  if (ci) return(.analytical_ci(res, n, are = .are_values[["iqr_scaled"]], level, "iqr_scaled"))
  res
}
