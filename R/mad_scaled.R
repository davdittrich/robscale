#' Scaled Median Absolute Deviation
#'
#' Computes the median absolute deviation from the median (or a user-supplied
#' center), scaled by a consistency constant for asymptotic normality under
#' the Gaussian model.
#'
#' @param x A numeric vector.
#' @param center Optional numeric scalar giving the central value from which to
#'   measure absolute deviations. Defaults to the median of \code{x}.
#' @param constant Consistency constant for asymptotic normality at the
#'   Gaussian. Defaults to \eqn{1/\Phi^{-1}(0.75) \approx 1.4826}{1/Phi^-1(0.75) ~ 1.4826}.
#'   Set to \code{1} for the raw median absolute deviation.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @details
#' The scaled MAD is defined as
#'
#' \deqn{\mathrm{MAD}_s(x) = C \cdot \mathrm{med}_i\, |x_i - \mathrm{med}(x)|}{MAD_s(x) = C * median(|x - median(x)|)}
#'
#' \strong{Statistical Properties.}
#' The MAD achieves a \bold{50\% breakdown point}, meaning it tolerates up to
#' half the data being contaminated. Its \bold{asymptotic relative efficiency
#' (ARE) is 36.8\%} compared to the sample standard deviation.
#'
#' \strong{Computational Performance.}
#' Unlike \code{\link[stats]{mad}}, this implementation uses O(n) selection
#' with adaptive algorithm dispatch: Floyd-Rivest for moderate n, pdqselect
#' for large n (threshold derived from per-core L2 cache size at startup),
#' and sorting networks for n \eqn{\le} 16.
#'
#' @return If \code{ci = FALSE} (default), a single numeric value: the scaled
#'   median absolute deviation. Returns \code{0} if \code{n = 1}; returns
#'   \code{NA} if \code{x} has length zero after removal of \code{NA}s.
#'   If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
#'
#' @seealso
#' \code{\link[stats]{mad}} for the base R implementation;
#' \code{\link{scale_robust}} for the unified dispatcher;
#' \code{\link{robScale}} for the M-estimate of scale that uses the MAD as
#' a starting value.
#'
#' @examples
#' mad_scaled(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' mad_scaled(x)                    # with consistency constant
#' mad_scaled(x, constant = 1)      # raw median absolute deviation
#'
#' # Supply a known center
#' mad_scaled(x, center = 4.0)
#'
#' # Asymptotic confidence interval
#' mad_scaled(x, ci = TRUE)
#'
#' @keywords univar robust
#' @export
mad_scaled <- function(x, center = NULL, constant = 1.482602218505602, na.rm = FALSE,
                       ci = FALSE, level = 0.95) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  n <- length(x)
  if (n == 0L) return(NA_real_)
  if (is.null(center)) {
    res <- .Call(`_robscale_mad_impl_auto`, x, constant)
  } else {
    res <- .Call(`_robscale_mad_impl_center`, x, center, constant)
  }
  if (ci) return(.analytical_ci(res, n, are = .are_values$mad, level, "mad_scaled"))
  res
}
