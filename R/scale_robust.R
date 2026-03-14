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
#' @param ci Logical. If \code{TRUE}, return confidence intervals alongside
#'   the point estimate. Single methods yield a \code{"robscale_ci"} object
#'   with analytical CIs; the ensemble yields a \code{"robscale_ensemble_ci"}
#'   object with both analytical and bootstrap CIs. Default: \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#' @param boot_method Bootstrap CI method for the ensemble. \code{"auto"}
#'   (default) selects BCa for n \eqn{\le} 200, percentile for
#'   n \eqn{\le} 5000, and parametric otherwise. Override with \code{"bca"},
#'   \code{"percentile"}, or \code{"parametric"}.
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
#' @return If \code{ci = FALSE} (default), a single numeric value: the scale
#'   estimate (\code{NA} when \code{n < 2}). If \code{ci = TRUE} with a single
#'   method, a \code{"robscale_ci"} object. If \code{ci = TRUE} with
#'   \code{method = "ensemble"}, a \code{"robscale_ensemble_ci"} object
#'   containing the ensemble estimate, bootstrap CI, and a table of
#'   per-estimator results.
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
#' # Ensemble with bootstrap CIs
#' scale_robust(x, ci = TRUE)
#'
#' # Single-method analytical CI
#' scale_robust(x, method = "sn", ci = TRUE)
#'
#' @keywords univar robust
#' @export
scale_robust <- function(x,
                         method = c("ensemble", "gmd", "sd", "mad", "iqr",
                                    "sn", "qn", "robScale"),
                         auto_switch = TRUE,
                         threshold = 20L,
                         n_boot = 200L,
                         na.rm = FALSE,
                         ci = FALSE,
                         level = 0.95,
                         boot_method = c("auto", "bca", "percentile",
                                         "parametric")) {
  method <- match.arg(method)
  boot_method <- match.arg(boot_method)
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  } else {
    if (anyNA(x)) {
      stop("There are NAs in the data yet na.rm is FALSE")
    }
  }
  n <- length(x)
  if (n < 2L) return(NA_real_)

  # Auto-switch to GMD for large n
  if (auto_switch && n >= threshold) {
    return(gmd(x, ci = ci, level = level))
  }

  # Non-ensemble methods
  if (method != "ensemble") {
    return(switch(method,
      gmd      = gmd(x, ci = ci, level = level),
      sd       = sd_c4(x, ci = ci, level = level),
      mad      = mad_scaled(x, ci = ci, level = level),
      iqr      = iqr_scaled(x, ci = ci, level = level),
      sn       = sn(x, ci = ci, level = level),
      qn       = qn(x, ci = ci, level = level),
      robScale = robScale(x, ci = ci, level = level)
    ))
  }

  # Ensemble without CI — fast path
  if (!ci) return(cpp_scale_ensemble(x, n_boot))

  # Ensemble with CI — determine bootstrap tier
  if (boot_method == "auto") {
    method_code <- if (n <= 200L) 0L else if (n <= 5000L) 1L else 2L
  } else {
    method_code <- switch(boot_method,
      bca = 0L, percentile = 1L, parametric = 2L
    )
  }

  raw <- cpp_scale_ensemble_ci(x, n_boot, level, method_code)

  # Analytical CIs for each component
  estimator_names <- c("sd_c4", "gmd", "mad_scaled", "iqr_scaled",
                       "sn", "qn", "robScale")
  are_values <- c(1.00, 0.98, 0.368, 0.37, 0.58, 0.82, 0.55)

  analytical_lower <- numeric(7L)
  analytical_upper <- numeric(7L)
  for (j in seq_len(7L)) {
    if (j == 1L) {
      ci_obj <- .chisq_ci(raw$estimates[j], n, level)
    } else {
      ci_obj <- .analytical_ci(raw$estimates[j], n, are_values[j],
                               level, estimator_names[j])
    }
    analytical_lower[j] <- ci_obj$ci[["lower"]]
    analytical_upper[j] <- ci_obj$ci[["upper"]]
  }

  boot_method_name <- c("bca", "percentile", "parametric")[method_code + 1L]

  estimators_df <- data.frame(
    estimator        = estimator_names,
    estimate         = raw$estimates,
    weight           = raw$weights,
    analytical_lower = analytical_lower,
    analytical_upper = analytical_upper,
    boot_lower       = raw$boot_lowers,
    boot_upper       = raw$boot_uppers,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      estimate    = raw$estimate,
      ci          = c(lower = raw$ci_lower, upper = raw$ci_upper),
      level       = level,
      method      = "ensemble",
      boot_method = boot_method_name,
      estimators  = estimators_df
    ),
    class = "robscale_ensemble_ci"
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
