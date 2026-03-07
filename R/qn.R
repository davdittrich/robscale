#' Robust Estimator of Scale Qn
#'
#' Computes the robust estimator of scale \eqn{Q_n} proposed by Rousseeuw and Croux (1993).
#'
#' @param x a numeric vector of observations.
#' @param constant consistency constant. Default is \eqn{2.2191}.
#' @param finite.corr logical; if \code{TRUE}, a finite-sample correction factor is applied.
#' @param na.rm logical; if \code{TRUE}, \code{NA} values are removed before computation.
#'
#' @return The \eqn{Q_n} estimator of scale.
#'
#' @details
#' \eqn{Q_n} is defined as \eqn{2.2191 \cdot \{ |x_i - x_j|; i < j \}_{(k)}} where \eqn{k \approx \binom{n}{2} / 4}.
#' It has a breakdown point of 50%, is more efficient than the MAD, and does not require a location estimate.
#'
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#'
#' @export
qn <- function(x, constant = 2.2191, finite.corr = TRUE, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  
  if (is.double(x)) {
    res <- C_qn_fast(x)
  } else if (is.integer(x)) {
    res <- C_qn_int_fast(x)
  } else {
    res <- C_qn_fast(as.double(x))
  }
  
  if (constant != 2.2191) {
    res <- res * (constant / 2.21914446598508)
  }
  
  if (!finite.corr) {
    n <- length(x)
    res <- res / robscale:::get_qn_factor_r(n)
  }
  
  res
}

# Helper for finite correction (internal)
get_qn_factor_r <- function(n) {
  if (n <= 100) {
    factors <- c(0.00000, 0.39954, 0.99386, 0.51333, 0.84412, 0.61224, 0.85886, 0.67000, 0.87359, 0.72007,
                0.88902, 0.75748, 0.90232, 0.78551, 0.91248, 0.80779, 0.92106, 0.82600, 0.92793, 0.84105,
                0.93380, 0.85367, 0.93894, 0.86441, 0.94303, 0.87372, 0.94680, 0.88186, 0.95009, 0.88901,
                0.95304, 0.89531, 0.95566, 0.90099, 0.95789, 0.90600, 0.96004, 0.91061, 0.96192, 0.91480,
                0.96361, 0.91852, 0.96522, 0.92200, 0.96668, 0.92515, 0.96802, 0.92809, 0.96923, 0.93085,
                0.97040, 0.93334, 0.97147, 0.93566, 0.97237, 0.93781, 0.97328, 0.93985, 0.97421, 0.94180,
                0.97496, 0.94355, 0.97573, 0.94525, 0.97648, 0.94687, 0.97710, 0.94837, 0.97773, 0.94978,
                0.97837, 0.95112, 0.97891, 0.95235, 0.97944, 0.95359, 0.97999, 0.95472, 0.98049, 0.95579,
                0.98090, 0.95677, 0.98138, 0.95781, 0.98179, 0.95871, 0.98216, 0.95967, 0.98255, 0.96051,
                0.98295, 0.96139, 0.98329, 0.96212, 0.98363, 0.96294, 0.98399, 0.96364, 0.98430, 0.96438)
    return(factors[n])
  }
  inv_n <- 1.0 / n
  if (n %% 2 == 1) {
    return(1.0 - 1.6022 * inv_n + 4.7453 * inv_n^2)
  } else {
    return(1.0 - 3.6741 * inv_n + 11.1030 * inv_n^2)
  }
}
