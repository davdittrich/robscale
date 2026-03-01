#' Robust M-Estimate of Scale
#'
#' Compute the robust M-estimate of scale for very small samples using the
#' \eqn{\rho}{rho} function of Rousseeuw & Verboven (2002).
#'
#' @param x A numeric vector.
#' @param loc Optional numeric scalar giving a known location.  When supplied,
#'   the observations are centered at \code{loc} and the minimum sample size
#'   for iteration is lowered from 4 to 3 (see \sQuote{Details}).
#' @param implbound Implosion bound: the smallest value the MAD is allowed to
#'   take before it is considered to have \dQuote{imploded} (collapsed to
#'   zero).  Defaults to \code{1e-4}.
#' @param na.rm Logical.  If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation.  If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param maxit Maximum number of multiplicative iterations.  Defaults to 80.
#' @param tol Convergence tolerance.  Iteration stops when the multiplicative
#'   update factor satisfies \eqn{|v - 1| \le \mathrm{tol}}{|v - 1| <= tol}.
#'   Defaults to \code{sqrt(.Machine$double.eps)}.
#'
#' @details
#' The scale estimator \eqn{S_n}{Sn} solves the M-estimating equation
#'
#' \deqn{\frac{1}{n}\sum_{i=1}^{n}\rho\!\left(\frac{x_i - T_n}{S_n}
#'   \right) = \beta}{mean(rho((x_i - Tn) / Sn)) = beta}
#'
#' where \eqn{T_n}{Tn} is fixed at the sample median, \eqn{\beta = 0.5}, and
#' \eqn{\rho} is a smooth rho function defined as the square of the logistic
#' psi (Rousseeuw & Verboven, 2002, Sec.\sspace{}4.2):
#'
#' \deqn{\rho_{\mathrm{log}}(x) = \psi_{\mathrm{log}}^2\!\left(
#'   \frac{x}{c}\right)}{rho(x) = psi(x / c)^2}
#'
#' with the tuning constant \eqn{c = 0.37394112142347236} chosen so that
#'
#' \deqn{\int\rho(u)\,d\Phi(u) = 0.5}{Int rho(u) dPhi(u) = 0.5}
#'
#' yielding a 50\% breakdown point.
#'
#' \strong{Iteration scheme.}
#' The equation is solved by multiplicative iteration (Rousseeuw & Verboven,
#' 2002, Eq.\sspace{}27):
#'
#' \deqn{S^{(k+1)} = S^{(k)} \cdot \sqrt{2 \cdot \frac{1}{n}\sum
#'   \psi_{\mathrm{log}}^2\!\left(\frac{x_i - T}{c \cdot
#'   S^{(k)}}\right)}}{S(k+1) = S(k) * sqrt(2 * mean(psi((x_i - T) /
#'   (c * S(k)))^2))}
#'
#' Starting value: \eqn{S^{(0)} = \mathrm{MAD}(x)}{S(0) = MAD(x)}.
#' The logistic psi values are computed via the algebraic identity
#' \eqn{\psi_{\mathrm{log}}(x) = \tanh(x/2)}{psi(x) = tanh(x/2)}.
#'
#' \strong{Decoupled estimation.}
#' Scale is estimated with location held fixed at
#' \eqn{\mathrm{med}(x)}{median(x)}, following the decoupled approach of
#' Rousseeuw & Verboven (2002, Sec.\sspace{}4.2).  This avoids the
#' positive-feedback instability of Huber's Proposal 2 in small samples.
#'
#' \strong{Known location.}
#' When \code{loc} is supplied, the observations are centered as
#' \eqn{x_i - \mu}{x_i - mu} and the initial scale is set to
#' \eqn{1.4826 \cdot \mathrm{med}(|x_i - \mu|)}{1.4826 * median(|x_i - mu|)}
#' rather than the MAD.  This lowers the minimum sample size from 4 to 3
#' (Rousseeuw & Verboven, 2002, Sec.\sspace{}5).
#'
#' \strong{Fallback.}
#' When \eqn{n} is below the minimum for iteration:
#' \itemize{
#'   \item if \eqn{\mathrm{MAD}(x) \le}{MAD(x) <=} \code{implbound}
#'     (implosion), the function returns \code{\link{adm}(x)};
#'   \item otherwise, it returns \eqn{\mathrm{MAD}(x)}{MAD(x)}.
#' }
#'
#' @return A single numeric value: the robust M-estimate of scale.
#'   Returns \code{NA} if \code{x} has length zero (after removal of
#'   \code{NA}s when \code{na.rm = TRUE}).
#'
#' @references
#' Rousseeuw, P. J. and Verboven, S. (2002) Robust estimation in very small
#' samples. \emph{Computational Statistics & Data Analysis}, \bold{40}(4),
#' 741--758. \doi{10.1016/S0167-9473(02)00078-6}
#'
#' @seealso
#' \code{\link{adm}} for the implosion fallback;
#' \code{\link[stats]{mad}} for the starting value and classical alternative;
#' \code{\link{robLoc}} for the companion location estimator.
#'
#' @examples
#' robScale(c(1:9))
#'
#' x <- c(1, 2, 3, 5, 7, 8)
#' robScale(x)
#' robScale(x, loc = 5)           # known location
#'
#' # Outlier resistance
#' x_clean <- c(2.0, 3.1, 2.7, 2.9, 3.3)
#' x_dirty <- replace(x_clean, 5, 100)
#' c(robScale(x_clean), robScale(x_dirty))   # barely moves
#' c(sd(x_clean), sd(x_dirty))               # destroyed
#'
#' # MAD implosion: identical values cause MAD = 0
#' robScale(c(5, 5, 5, 5, 6))     # falls back to adm()
#'
#' @keywords univar robust
#' @export
robScale <- function(x, loc = NULL, implbound = 1e-4, na.rm = FALSE,
                     maxit = 80L, tol = sqrt(.Machine$double.eps)) {
  if (na.rm) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  if (length(x) == 0L) return(NA_real_)
  has_loc <- !is.null(loc)
  loc_val <- if (has_loc) loc else 0.0
  rob_scale_impl(x, has_loc, loc_val, implbound, maxit, tol)
}
