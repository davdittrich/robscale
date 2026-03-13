#' Robust Ensemble Scale Estimation
#'
#' Unified dispatcher for robust scale estimation. Automatically selects
#' between a variance-weighted ensemble of 7 estimators (for small samples)
#' and the Gini Mean Difference (for large samples), or returns a specific
#' estimator by name.
#'
#' @param x A numeric vector of observations.
#' @param method Character string specifying the estimation method.
#'   Options: \code{"ensemble"} (default), \code{"gmd"}, \code{"sd"},
#'   \code{"mad"}, \code{"iqr"}, \code{"sn"}, \code{"qn"}, \code{"robScale"}.
#' @param auto_switch Logical. If \code{TRUE} (default), automatically uses
#'   GMD for \code{n >= threshold}.
#' @param threshold Integer. Sample size at which the switch to GMD occurs.
#'   Default: 20 (research-backed; see \sQuote{Details}).
#' @param n_boot Integer. Number of bootstrap replicates for the ensemble
#'   weighting. Default: 200.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. Default: \code{FALSE}.
#'
#' @details
#' \strong{Ensemble method.}
#' When \code{method = "ensemble"}, the function computes a
#' variance-weighted combination of 7 scale estimators:
#' \enumerate{
#'   \item \code{\link{sd_c4}} --- unbiased standard deviation
#'   \item \code{\link{gmd}} --- Gini mean difference
#'   \item \code{\link{mad_scaled}} --- median absolute deviation
#'   \item \code{\link{iqr_scaled}} --- scaled interquartile range
#'   \item \code{\link{sn}} --- Sn estimator of Rousseeuw & Croux
#'   \item \code{\link{qn}} --- Qn estimator of Rousseeuw & Croux
#'   \item \code{\link{robScale}} --- logistic M-estimator of scale
#' }
#' The weights are determined by bootstrap resampling: each estimator's
#' inverse variance across \code{n_boot} resamples determines its
#' contribution. Estimators with lower sampling variance receive higher
#' weight.
#'
#' \strong{Automatic switching.}
#' When \code{auto_switch = TRUE} and \code{n >= threshold}, the function
#' returns \code{\link{gmd}(x)} directly. The GMD achieves 98\% asymptotic
#' relative efficiency at the Gaussian while being computationally cheaper
#' than the ensemble, making it the optimal choice for moderate-to-large
#' samples.
#'
#' \strong{Individual methods.}
#' When a specific method is requested, \code{scale_robust} dispatches to the
#' corresponding function with default parameters.
#'
#' @return A single numeric value: the scale estimate.
#'   Returns \code{NA} if \code{n < 2}.
#'
#' @seealso
#' Individual estimators: \code{\link{sd_c4}}, \code{\link{gmd}},
#' \code{\link{mad_scaled}}, \code{\link{iqr_scaled}}, \code{\link{sn}},
#' \code{\link{qn}}, \code{\link{robScale}};
#' \code{\link{get_consistency_constant}} for the underlying constants.
#'
#' @examples
#' x <- c(1, 2, 3, 5, 7, 8)
#' scale_robust(x)                          # ensemble (n < 20)
#' scale_robust(x, method = "qn")           # specific method
#'
#' set.seed(42)
#' y <- rnorm(50)
#' scale_robust(y)                          # auto-switches to GMD (n >= 20)
#' scale_robust(y, auto_switch = FALSE)     # forces ensemble
#'
#' @keywords univar robust
#' @export
scale_robust <- function(x,
                         method = c("ensemble", "gmd", "sd", "mad", "iqr",
                                    "sn", "qn", "robScale"),
                         auto_switch = TRUE,
                         threshold = 20L,
                         n_boot = 200L,
                         na.rm = FALSE) {
  method <- match.arg(method)
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  n <- length(x)
  if (n < 2L) return(NA_real_)

  if (auto_switch && n >= threshold) return(gmd(x))

  switch(method,
    ensemble = cpp_scale_ensemble(x, n_boot),
    gmd      = gmd(x),
    sd       = sd_c4(x),
    mad      = mad_scaled(x),
    iqr      = iqr_scaled(x),
    sn       = sn(x),
    qn       = qn(x),
    robScale = robScale(x)
  )
}

#' Get Consistency Constant
#'
#' Returns the consistency constant or finite-sample correction factor for a
#' given scale estimator and sample size.
#'
#' @param method Character string specifying the estimator.
#'   Options: \code{"c4"}, \code{"gmd"}, \code{"mad"}, \code{"iqr"},
#'   \code{"sn"}, \code{"qn"}.
#' @param n Integer. The sample size (used for \code{"c4"}, \code{"sn"},
#'   and \code{"qn"} which have sample-size-dependent factors).
#'
#' @return A single numeric value: the consistency constant or correction
#'   factor.
#'
#' @details
#' For \code{"c4"}, \code{"sn"}, and \code{"qn"}, the returned value depends
#' on \code{n}. For \code{"gmd"}, \code{"mad"}, and \code{"iqr"}, the
#' returned value is the asymptotic constant (independent of \code{n}).
#'
#' @examples
#' get_consistency_constant("mad")      # 1.4826
#' get_consistency_constant("c4", 5)    # c4(5)
#' get_consistency_constant("qn", 10)   # finite-sample Qn factor
#'
#' @keywords univar
#' @export
get_consistency_constant <- function(method, n = NULL) {
  method <- match.arg(method, c("c4", "gmd", "mad", "iqr", "sn", "qn"))
  switch(method,
    c4  = {
      if (is.null(n)) stop("n is required for c4")
      # Use the same formula as the C++ c4_factor
      exp(0.5 * log(2 / (n - 1)) + lgamma(n / 2) - lgamma((n - 1) / 2))
    },
    gmd = 0.886226925452758,
    mad = 1.482602218505602,
    iqr = 0.741301109252801,
    sn  = {
      if (is.null(n)) stop("n is required for sn finite-sample correction")
      C_get_sn_factor(as.integer(n))
    },
    qn  = {
      if (is.null(n)) stop("n is required for qn finite-sample correction")
      C_get_qn_factor(as.integer(n))
    }
  )
}
