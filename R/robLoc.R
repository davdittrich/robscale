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
#' \bold{asymptotic relative efficiency (ARE) of 0.98} compared to the sample
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
#' The denominator of the Newton step is the \emph{observed Fisher information}
#'
#' \deqn{\hat{I}(t) = \sum_{i=1}^{n}
#'   \operatorname{sech}^2\!\left(\frac{x_i - t}{2S_n}\right),}{
#'   I.hat(t) = sum(sech^2((x_i - t) / (2 * Sn))),}
#'
#' recomputed at each iteration using the current estimate \eqn{t}{t}.  This
#' distinguishes the implementation from IRLS (iteratively reweighted least
#' squares), which uses a fixed expected-information denominator and converges
#' linearly.  Using the observed Hessian yields true Newton--Raphson: the
#' fixed-point derivative at the solution satisfies \eqn{T'(t^*) = 0}{T'(t*)
#' = 0}, giving \emph{quadratic} local convergence.  On typical data,
#' convergence is achieved in two to four iterations; the \code{maxit} bound
#' of 80 is a conservative safety limit.
#'
#' \strong{Performance and SIMD.}
#' The C++ kernel dispatches \code{tanh} to the fastest available backend:
#' Apple Accelerate on macOS, glibc libmvec (AVX-512 8-wide or AVX2 4-wide)
#' on Linux x86-64, SLEEF when libmvec is absent, or \code{#pragma omp simd}
#' as a portable fallback. A fused AVX2 kernel accumulates \eqn{\psi_i} and
#' \eqn{\mathrm{d}\psi_i} in a single pass, halving memory reads.
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
#' \code{\link{robScale}} for the companion scale estimator;
#' \code{\link{qn}} and \code{\link{sn}} for high-efficiency scale
#' estimators.
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
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (!is.null(scale)) {
    if (!is.numeric(scale) || length(scale) != 1L)
      stop("'scale' must be a single numeric value")
  }
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  has_scale <- !is.null(scale)
  scale_val <- if (has_scale) scale else 0.0
  .Call(`_robscale_rob_loc_impl`, x, has_scale, scale_val, maxit, tol)
}
