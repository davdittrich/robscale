#' Robust M-Estimate of Scale
#'
#' Computes the robust M-estimate of scale for very small samples using the
#' \eqn{\rho}{rho} function of Rousseeuw & Verboven (2002).
#'
#' @param x A numeric vector.
#' @param loc Optional numeric scalar giving a known location. When supplied,
#'   the observations are centered at \code{loc} and the minimum sample size
#'   for iteration is lowered from 4 to 3 (see \sQuote{Details}).
#' @param fallback Character string specifying the fallback behavior when the MAD
#'   collapses to zero or the sample size is too small for iteration.
#'   Must be one of \code{"adm"} (default) or \code{"na"}. See \sQuote{Details}.
#' @param implbound Numeric scalar specifying the threshold for MAD implosion.
#'   Defaults to \code{1e-4}. Passing a value of 0 disables implosion checks.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. If \code{FALSE} (the default), the presence
#'   of any \code{NA} raises an error.
#' @param maxit Maximum number of Newton--Raphson iterations.
#'   Defaults to 80.
#' @param tol Convergence tolerance. Iteration stops when the relative
#'   change in the scale estimate falls below \code{tol}. Defaults to
#'   \code{sqrt(.Machine$double.eps)}.
#' @param ci Logical. If \code{TRUE}, return a \code{"robscale_ci"} object
#'   with the point estimate and asymptotic confidence interval.
#'   Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#'
#' @details
#' The M-estimate of scale \eqn{S_n}{Sn} is defined as the solution to the
#' estimating equation
#'
#' \deqn{\frac{1}{n} \sum_{i=1}^{n} \rho \left( \frac{x_i - T_n}{S_n} \right) = \beta}{mean(rho((x_i - Tn) / Sn)) = beta}
#'
#' where the location \eqn{T_n}{Tn} is fixed at the sample median, \eqn{\beta = 0.5}
#' is the expected value of \eqn{\rho} under the Gaussian model, and \eqn{\rho}
#' is a smooth rho function (Rousseeuw & Verboven, 2002, Sec.\sspace{}4.2):
#'
#' \deqn{\rho_{\mathrm{log}}(x) = \psi_{\mathrm{log}}^2 \left( \frac{x}{c} \right)}{rho(x) = psi(x / c)^2}
#'
#' The tuning constant \eqn{c = 0.373941121} is chosen to satisfy
#' \eqn{E_\Phi[\rho(u)] = 0.5}.
#'
#' \strong{Statistical Properties.}
#' This estimator is designed for high robustness and efficiency. It achieves a
#' \bold{50\% breakdown point}, meaning the estimate remains reliable even if
#' half the sample is contaminated by outliers. At the Gaussian distribution,
#' the logistic M-estimator of scale achieves an \bold{asymptotic relative
#' efficiency (ARE) of 0.55} compared to the sample standard deviation.
#'
#' \strong{Numerical Computation.}
#' The M-scale estimating equation \eqn{n^{-1}\sum\rho(u_i) = 1/2} is solved
#' by Newton--Raphson iteration, starting from the MAD. Each step computes
#' \eqn{u_i = (x_i - T)/(2cS)} and accumulates two sums in a single pass:
#' \eqn{\sum \tanh^2(u_i)}{sum(tanh^2(u_i))} (the rho sum) and
#' \eqn{\sum u_i \tanh(u_i)\operatorname{sech}^2(u_i)}{sum(u_i*tanh(u_i)*sech^2(u_i))}
#' (the derivative sum). The NR update is
#' \deqn{\Delta S = S \cdot \frac{\bar\rho - 1/2}{(2/n)\sum u_i \tanh(u_i)\operatorname{sech}^2(u_i)}}{dS = S * (mean_rho - 0.5) / ((2/n) * sum(u*tanh(u)*sech^2(u)))}
#' with convergence test \eqn{|\Delta S|/S \le \mathrm{tol}}{|dS|/S <= tol}.
#' When the derivative sum degenerates, the iteration falls back to a
#' multiplicative half-step. Quadratic convergence yields 3--4 iterations
#' on typical data. Location is held fixed at the sample median, following
#' the decoupled approach of Rousseeuw and Verboven (2002) that avoids the
#' positive-feedback instabilities of simultaneous location--scale estimation.
#'
#' \strong{Performance and SIMD.}
#' The C++ kernel dispatches \code{tanh} to the fastest available backend:
#' Apple Accelerate on macOS, glibc libmvec (AVX-512 8-wide or AVX2 4-wide)
#' on Linux x86\_64, SLEEF when libmvec is absent, or \code{#pragma omp simd}
#' as a portable fallback.
#'
#' \strong{Known location.}
#' When \code{loc} is supplied, the observations are centered as
#' \eqn{x_i - \mu}{x_i - mu} and the initial scale is set to
#' \eqn{1.4826 \cdot \mathrm{med}(|x_i - \mu|)}{1.4826 * median(|x_i - mu|)}
#' rather than the MAD.  This lowers the minimum sample size from 4 to 3
#' (Rousseeuw & Verboven, 2002, Sec.\sspace{}5).
#'
#' \strong{Fallback Mechanism and Implosion.}
#' Robust scale estimators like the MAD can "implode" (collapse to zero) if more
#' than 50% of the sample observations are identical. When the MAD collapses
#' or the sample size is too small for reliable iteration (\eqn{n < 4}, or
#' \eqn{n < 3} if location is known):
#' \itemize{
#'   \item If \code{fallback = "adm"} (default), the function returns the
#'     scaled Average Distance to the Median (\code{\link{adm}}). The ADM
#'     is highly resistant to implosion (breakdown point \eqn{(n-1)/n}).
#'   \item If \code{fallback = "na"}, the function returns \code{NA}.
#' }
#'
#' @return If \code{ci = FALSE} (default), a single numeric value: the robust
#'   M-estimate of scale. Returns \code{NA} if \code{x} has length zero (after
#'   removal of \code{NA}s when \code{na.rm = TRUE}) or if the MAD collapses
#'   and \code{fallback = "na"}.
#'   If \code{ci = TRUE}, an object of class \code{"robscale_ci"}.
#'
#' @references
#' Rousseeuw, P. J. and Verboven, S. (2002) Robust estimation in very small
#' samples. \emph{Computational Statistics & Data Analysis}, \bold{40}(4),
#' 741--758. \doi{10.1016/S0167-9473(02)00078-6}
#'
#' @seealso
#' \code{\link{adm}} for the implosion fallback;
#' \code{\link[stats]{mad}} for the starting value and classical alternative;
#' \code{\link{robLoc}} for the companion location estimator;
#' \code{\link{qn}} and \code{\link{sn}} for high-efficiency scale
#' estimators.
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
#' # Asymptotic confidence interval
#' robScale(x, ci = TRUE)
#'
#' @keywords univar robust
#' @export
robScale <- function(x, loc = NULL, fallback = c("adm", "na"),
                     implbound = 1e-4, na.rm = FALSE,
                     maxit = 80L, tol = sqrt(.Machine$double.eps),
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

  fallback <- match.arg(fallback)
  fallback_code <- if (fallback == "adm") 0L else 1L

  if (!is.null(loc)) {
    if (!is.numeric(loc) || length(loc) != 1L)
      stop("'loc' must be a single numeric value")
  }

  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }

  has_loc <- !is.null(loc)
  loc_val <- if (has_loc) loc else 0.0
  res <- .Call(`_robscale_rob_scale_impl`, x, has_loc, loc_val, implbound, maxit, tol, fallback_code)
  if (ci) return(.analytical_ci(res, n, are = .are_values[["robScale"]], level, "robScale"))
  res
}
