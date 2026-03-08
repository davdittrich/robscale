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
#'
#' @details
#' The average distance to the median (ADM) is defined as
#'
#' \deqn{\mathrm{ADM}(x) \;=\; C \cdot \frac{1}{n}\sum_{i=1}^{n}
#'   |x_i - \mathrm{med}(x)|}{ADM(x) = C * mean(|x - med(x)|)}
#'
#' where \eqn{C} is the consistency constant and \eqn{\mathrm{med}(x)}{med(x)}
#' is the sample median.  When \code{center} is supplied it replaces the
#' median.
#'
#' The default constant \eqn{C = \sqrt{\pi/2}}{C = sqrt(pi/2)} makes the ADM a
#' consistent estimator of the standard deviation under the normal distribution.
#' In large samples the ADM converges to \eqn{\sigma} at the Gaussian;
#' however, this asymptotic property may \strong{not} hold at the very small
#' sample sizes (\eqn{n = 3}--\eqn{8}) for which this package is primarily
#' intended.
#'
#' The ADM is \emph{not} robust against outliers in the explosion sense:
#' its explosion breakdown point is \eqn{1/n}.  However, it is highly
#' resistant to implosion, with an implosion breakdown point of
#' \eqn{(n-1)/n}{(n - 1) / n}.  It therefore serves as the fallback scale
#' estimator in \code{\link{robScale}} when the MAD collapses to zero.
#'
#' \strong{Performance.}
#' This implementation uses a C++17 \eqn{O(n)} selection strategy (sorting networks
#' for \eqn{n \le 8}, introselect for larger \eqn{n}), typically yielding a 10x
#' speedup over pure-R implementations like \code{revss}.
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
adm <- function(x, center, constant = 1.2533141373155001, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (length(x) == 0L) return(NA_real_)
  if (missing(center)) {
    adm_impl_auto(x, constant)
  } else {
    adm_impl(x, center, constant)
  }
}
