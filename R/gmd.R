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
#' alternative in this package. Its breakdown point is approximately 29.3\%.
#'
#' @return A single numeric value: the scaled Gini mean difference.
#'   Returns \code{0} if \code{n < 2}; returns \code{NA} if \code{x} has
#'   length zero after removal of \code{NA}s.
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
#' @keywords univar robust
#' @export
gmd <- function(x, constant = 0.886226925452758, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (length(x) == 0L) return(NA_real_)
  gmd_impl(x, constant)
}
