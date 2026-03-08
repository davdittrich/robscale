#' Robust M-Estimate of Location
#'
#' Computes the robust M-estimate of location for very small samples using the
#' logistic \eqn{\psi}{psi} function of Rousseeuw & Verboven (2002).
#'
#' @param x A numeric vector.
#' @param scale Optional numeric scalar giving a known scale.  When supplied,
#'   the MAD is replaced by this value and the minimum sample size for
#'   iteration is lowered from 4 to 3 (see \sQuote{Details}).
#' @param na.rm Logical.  If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation.  If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param maxit Maximum number of Newton--Raphson iterations.  Defaults to 80.
#' @param tol Convergence tolerance.  Iteration stops when the absolute
#'   Newton step falls below \code{tol}.  Defaults to
#'   \code{sqrt(.Machine$double.eps)}.
#'
#' @details
#' The M-estimate of location \eqn{T_n}{Tn} is defined as the solution to the
#' estimating equation
#'
#' \deqn{\sum_{i=1}^{n}\psi_{\mathrm{log}}
#'   \!\left(\frac{x_i - T_n}{S_n}\right) = 0}{
#'   sum(psi((x_i - Tn) / Sn)) = 0}
#'
#' where \eqn{S_n}{Sn} is a fixed auxiliary scale (defaulting to the MAD) and
#' \eqn{\psi_{\mathrm{log}}}{psi} is the logistic psi function:
#'
#' \deqn{\psi_{\mathrm{log}}(x) = \frac{e^x - 1}{e^x + 1}
#'   = \tanh(x/2)}{psi(x) = (exp(x) - 1) / (exp(x) + 1) = tanh(x / 2)}
#'
#' \strong{Statistical Properties.}
#' The logistic psi function is bounded, smooth (\eqn{C^\infty}{C-inf}), and
#' strictly monotone. These properties ensure that the resulting M-estimator is
#' both robust to outliers and numerically stable. At the Gaussian distribution,
#' the logistic M-estimator of location achieves high efficiency, with an
#' \bold{asymptotic relative efficiency (ARE) of 0.95} compared to the sample
#' mean.
#'
#' \strong{Small-Sample Strategy.}
#' Following Rousseeuw & Verboven (2002), location and scale are estimated
#' separately. In \code{robLoc}, the auxiliary scale \eqn{S_n}{Sn} remains fixed
#' throughout the Newton--Raphson iteration. This "decoupled" approach avoids
#' the instabilities often encountered in small samples when using simultaneous
#' location--scale iteration (e.g., Huber's Proposal 2).
#'
#' \strong{Numerical Computation.}
#' The estimating equation is solved via Newton--Raphson iteration starting from
#' the sample median. Because the derivative of the logistic psi satisfies
#' \eqn{\psi'(x) = \frac{1}{2}(1 - \psi^2(x))}{\psi'(x) = 1/2 * (1 - \psi(x)^2)},
#' the Newton step is computationally efficient, requiring no additional
#' transcendental calls beyond the \code{tanh} evaluations used for the psi
#' function itself.
#'
#' \strong{Performance and SIMD.}
#' The underlying C++ core utilizes platform-specific SIMD backends (SLEEF on
#' Linux, Apple Accelerate on macOS) to vectorize the \code{tanh} evaluations.
#' This architectural choice delivers substantial performance gains,
#' particularly for large-scale or high-throughput workflows.
#'
#' \strong{Fallback Mechanism.}
#' For extremely small samples where iteration may be unreliable, the function
#' returns the \code{\link{median}} directly:
#' \itemize{
#'   \item If scale is unknown: \eqn{n < 4} (since the MAD of 3 points is
#'     insufficiently robust).
#'   \item If scale is supplied: \eqn{n < 3}.
#' }
#'
#' @return A single numeric value: the robust M-estimate of location.
#'   Returns \code{NA} if \code{x} has length zero (after removal of
#'   \code{NA}s when \code{na.rm = TRUE}).
#'
#' @references
#' Rousseeuw, P. J. and Verboven, S. (2002) Robust estimation in very small
#' samples. \emph{Computational Statistics & Data Analysis}, \bold{40}(4),
#' 741--758. \doi{10.1016/S0167-9473(02)00078-6}
#'
#' @seealso
#' \code{\link{median}} for the starting value;
#' \code{\link[stats]{mad}} for the auxiliary scale;
#' \code{\link{robScale}} for the companion scale estimator.
#'
#' @examples
#' robLoc(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' robLoc(x)
#'
#' # Known scale lowers the minimum sample size to 3
#' robLoc(c(1, 2, 3), scale = 1.5)
#'
#' # Outlier resistance
#' x_clean <- c(2.0, 3.1, 2.7, 2.9, 3.3)
#' x_dirty <- replace(x_clean, 5, 100)
#' c(robLoc(x_clean), robLoc(x_dirty))   # barely moves
#' c(mean(x_clean), mean(x_dirty))       # destroyed
#'
#' @keywords univar robust
#' @export
robLoc <- function(x, scale = NULL, na.rm = FALSE, maxit = 80L,
                   tol = sqrt(.Machine$double.eps)) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (length(x) == 0L) return(NA_real_)
  has_scale <- !is.null(scale)
  scale_val <- if (has_scale) scale else 0.0
  rob_loc_impl(x, has_scale, scale_val, maxit, tol)
}
