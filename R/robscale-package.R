#' Faster Robustness: SIMD-Accelerated Estimation of Location and Scale
#'
#' A comprehensive suite of high-performance robust scale and location
#' estimators. The package provides:
#'
#' \strong{Scale estimators:}
#' \itemize{
#'   \item \code{\link{qn}}, \code{\link{sn}} --- Rousseeuw & Croux (1993) with
#'     SIMD and TBB parallelism
#'   \item \code{\link{robScale}} --- logistic M-estimator with Newton-Raphson
#'     and SIMD-accelerated tanh
#'   \item \code{\link{gmd}} --- Gini mean difference (98\% ARE)
#'   \item \code{\link{mad_scaled}}, \code{\link{iqr_scaled}} --- scaled MAD and
#'     IQR with O(n) selection
#'   \item \code{\link{sd_c4}} --- unbiased standard deviation with c4 correction
#'   \item \code{\link{adm}} --- average distance to the median (implosion fallback)
#' }
#'
#' \strong{Location estimator:}
#' \itemize{
#'   \item \code{\link{robLoc}} --- logistic M-estimator of location
#' }
#'
#' \strong{Unified API:}
#' \itemize{
#'   \item \code{\link{scale_robust}} --- adaptive dispatcher with
#'     variance-weighted ensemble for small samples, automatic GMD switching
#'     for larger samples
#' }
#'
#' @return None, this is a package documentation object.
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
