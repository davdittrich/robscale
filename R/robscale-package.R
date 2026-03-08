#' Faster Robustness: SIMD-Accelerated Estimation of Location and Scale
#'
#' Robust estimation ensures statistical reliability in data 
#' contaminated by outliers. Yet, computational bottlenecks in existing R 
#' implementations frequently obstruct both very small sample analysis and 
#' large-scale processing. 'robscale' resolves these inefficiencies by 
#' providing high-performance C++17 implementations of logistic M-estimators 
#' and the Qn and Sn scale estimators. By leveraging platform-specific SIMD 
#' vectorization and TBB parallelism, the package delivers speedups of 
#' 11–39x for small samples and up to 10x for massive datasets. These 
#' performance gains enable the integration of robust statistics into modern, 
#' time-critical computational workflows.
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
