#' Fast Robust Estimation of Location and Scale
#'
#' Fast C++17 implementation of the robust location and scale estimators
#' of Rousseeuw & Verboven (2002) for very small samples, and the Qn and Sn
#' estimators of Rousseeuw and Croux (1993) for general samples.
#' Exploits the tanh identity for the logistic psi function and provides
#' platform-vectorized transcendental evaluation (SLEEF, Accelerate, OpenMP SIMD).
#' Parallized via RcppParallel/TBB for large samples.
#'
#' @keywords internal
#' @references
#' Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the Median Absolute Deviation.
#' \emph{Journal of the American Statistical Association}, \bold{88}, 1273--1283.
#'
#' Rousseeuw, P. J. and Verboven, S. (2002) Robust estimation in very small
#' samples. \emph{Computational Statistics & Data Analysis}, \bold{40}(4),
#' 741--758. \doi{10.1016/S0167-9473(02)00078-6}
#'
#' Nair, K. R. (1947) A Note on the Mean Deviation from the Median.
#' \emph{Biometrika}, \bold{34}(3/4), 360--362. \doi{10.2307/2332448}
"_PACKAGE"

#' @useDynLib robscale, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL
