#' Robust M-Estimate of Location
#'
#' Compute the robust M-estimate of location for very small samples using the
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
#' The location estimator \eqn{T_n}{Tn} is the solution to the
#' M-estimating equation
#'
#' \deqn{\frac{1}{n}\sum_{i=1}^{n}\psi_{\mathrm{log}}
#'   \!\left(\frac{x_i - T_n}{S_n}\right) = 0}{
#'   mean(psi((x_i - Tn) / Sn)) = 0}
#'
#' where \eqn{S_n}{Sn} is a fixed auxiliary scale estimate and
#' \eqn{\psi_{\mathrm{log}}}{psi} is the logistic psi function (Rousseeuw &
#' Verboven 2002, Eq.\sspace{}23):
#'
#' \deqn{\psi_{\mathrm{log}}(x) = \frac{e^x - 1}{e^x + 1}
#'   = \tanh(x/2)}{psi(x) = (exp(x) - 1) / (exp(x) + 1) = tanh(x / 2)}
#'
#' This function is bounded in \eqn{(-1, 1)}, smooth (\eqn{C^\infty}{C-inf}),
#' and strictly monotone.  Boundedness provides robustness against outliers;
#' smoothness avoids the corner artifacts of Huber's \eqn{\psi}{psi} at small
#' \eqn{n}.
#'
#' \strong{Iteration scheme.}
#' The estimating equation is solved by Newton--Raphson iteration.  The
#' derivative of the logistic psi satisfies \eqn{\psi'(x) =
#' 1 - \psi^2(x)}{psi'(x) = 1 - psi(x)^2}, so the Newton step requires
#' no additional transcendental function evaluations beyond those already
#' computed for the numerator.  Starting value:
#' \eqn{T^{(0)} = \mathrm{med}(x)}{T(0) = median(x)}.  Auxiliary scale:
#' \eqn{S = \mathrm{MAD}(x)}{S = MAD(x)} unless \code{scale} is supplied.
#'
#' \strong{Performance and SIMD.}
#' This C++17 implementation uses platform-specific SIMD vectorization to accelerate
#' transcendental evaluations (\code{tanh}). On Linux, it uses the \bold{SLEEF}
#' library (targeting AVX2, AVX512, or NEON); on macOS, it uses \bold{Apple Accelerate}.
#' These optimizations reduce execution time by a factor of 15--30 compared to
#' interpreted R implementations like \code{revss}.
#'
#' \strong{Decoupled estimation.}
#' Location and scale are estimated separately: \code{robLoc} holds the
#' auxiliary scale fixed at \eqn{\mathrm{MAD}(x)}{MAD(x)} throughout
#' iteration, following the decoupled approach of Rousseeuw & Verboven (2002,
#' Sec.\sspace{}4.1).  This avoids the positive-feedback instability of
#' simultaneous location--scale iteration (Huber's Proposal 2) in small
#' samples.
#'
#' \strong{Fallback.}
#' When the sample is too small for reliable iteration the function returns
#' the \code{median} directly:
#' \itemize{
#'   \item \eqn{n < 4} when \code{scale} is unknown (the MAD is unreliable
#'     at \eqn{n = 3});
#'   \item \eqn{n < 3} when \code{scale} is known.
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
