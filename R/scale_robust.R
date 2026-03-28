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
#'   GMD when \code{method = "ensemble"} and \code{n >= threshold}. Has no
#'   effect on named methods: \code{method = "qn"} always returns the Qn
#'   estimator regardless of sample size.
#' @param threshold Integer. Sample size at which the ensemble auto-switches
#'   to GMD. Default: 20 (research-backed; see \sQuote{Details}).
#' @param n_boot Integer. Number of bootstrap replicates for the ensemble
#'   weighting or bootstrap CI. Default: 200.
#' @param na.rm Logical. If \code{TRUE}, \code{NA} values are stripped from
#'   \code{x} before computation. Default: \code{FALSE}.
#' @param ci Logical. If \code{TRUE}, return confidence intervals alongside
#'   the point estimate. Single methods yield a \code{"robscale_ci"} object;
#'   the ensemble yields a \code{"robscale_ensemble_ci"} object. Default:
#'   \code{FALSE}.
#' @param level Confidence level for the interval (default 0.95).
#' @param boot_method CI method. For single named methods: \code{"auto"} and
#'   \code{"analytical"} (default) return the closed-form ARE-based interval
#'   (\code{sd} uses chi-squared; all others use the normal approximation).
#'   \code{"bca"}, \code{"percentile"}, and \code{"parametric"} return a
#'   bootstrap CI via \code{n_boot} resamples. For the ensemble: \code{"auto"}
#'   selects BCa for n \eqn{\le} 200, percentile for n \eqn{\le} 5000, and
#'   parametric otherwise. \code{"analytical"} is not supported for
#'   \code{method = "ensemble"}.
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
#' When \code{auto_switch = TRUE} and \code{method = "ensemble"} and
#' \code{n >= threshold}, the function returns \code{\link{gmd}(x)} directly.
#' The GMD achieves 98\% asymptotic relative efficiency at the Gaussian while
#' being computationally cheaper than the ensemble. Named methods (e.g.\
#' \code{method = "qn"}) are always dispatched as requested; \code{auto_switch}
#' never overrides an explicit method choice.
#'
#' \strong{Individual methods.}
#' When a specific method is requested, \code{scale_robust} bypasses the
#' ensemble and calls the corresponding C++ entry point directly. With
#' \code{ci = TRUE}, the default \code{boot_method = "auto"} returns the
#' analytical interval: chi-squared for \code{"sd"}, ARE-based normal
#' approximation for all others. Pass \code{boot_method = "bca"},
#' \code{"percentile"}, or \code{"parametric"} to obtain a bootstrap CI
#' instead.
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
#' scale_robust(x, method = "qn")           # specific method (always qn)
#'
#' set.seed(42)
#' y <- rnorm(50)
#' scale_robust(y)                          # ensemble auto-switches to GMD
#' scale_robust(y, method = "qn")           # qn regardless of n
#' scale_robust(y, auto_switch = FALSE)     # forces ensemble
#'
#' # Analytical CI for a named method (default)
#' scale_robust(x, method = "sn", ci = TRUE)
#'
#' # Bootstrap CI for a named method
#' scale_robust(x, method = "qn", ci = TRUE, boot_method = "percentile")
#'
#' # Ensemble with bootstrap CIs
#' scale_robust(x, ci = TRUE)
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
                                         "parametric", "analytical")) {
  method <- match.arg(method)

  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)

  boot_method <- match.arg(boot_method)

  # Resolve effective estimator: auto_switch only applies to the ensemble path.
  # Named methods (gmd, sd, …) are always dispatched as requested.
  effective_method <- method
  if (method == "ensemble" && auto_switch && n >= threshold)
    effective_method <- "gmd"

  if (ci) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
      stop("'level' must be a single numeric value in (0, 1)")
  }

  # "analytical" CI is not defined for the ensemble (which has no single ARE).
  if (boot_method == "analytical" && effective_method == "ensemble")
    stop('"analytical" CI is not available for method = "ensemble"')

  # ── Individual estimator path ────────────────────────────────────────
  if (effective_method != "ensemble") {
    est <- switch(effective_method,
      gmd      = gmd_impl(x, 0.886226925452758),
      sd       = sd_c4_impl(x),
      mad      = mad_impl_auto(x, 1.482602218505602),
      iqr      = iqr_impl(x, 0.741301109252801),
      sn       = C_sn_fast(x),
      qn       = C_qn_fast(x),
      robScale = C_rob_scale_fast(x)
    )

    if (!ci) return(est)

    # Canonical name used by .are_values and print output.
    result_method <- switch(effective_method,
      sd  = "sd_c4",
      mad = "mad_scaled",
      iqr = "iqr_scaled",
      effective_method       # gmd, sn, qn, robScale pass through unchanged
    )

    if (boot_method %in% c("auto", "analytical")) {
      ci_obj <- if (effective_method == "sd") {
        .chisq_ci(est, n, level)
      } else {
        .analytical_ci(est, n, .are_values[[result_method]], level,
                       result_method)
      }
      ci_obj$boot_method <- "analytical"
      return(ci_obj)
    }

    # Bootstrap CI for individual estimator
    estimator_id <- match(effective_method,
                          c("gmd", "sd", "mad", "iqr", "sn", "qn",
                            "robScale")) - 1L
    method_code  <- switch(boot_method, bca = 0L, percentile = 1L,
                            parametric = 2L)
    bounds <- cpp_single_estimator_ci_bounds(x, est, estimator_id, n_boot,
                                             level, method_code)
    return(structure(
      list(
        estimate    = est,
        ci          = c(lower = bounds$ci_lower, upper = bounds$ci_upper),
        level       = level,
        method      = result_method,
        boot_method = boot_method
      ),
      class = "robscale_ci"
    ))
  }

  # ── Ensemble path ────────────────────────────────────────────────────
  if (!ci) return(cpp_scale_ensemble(x, n_boot))

  # Determine bootstrap tier for ensemble CI
  if (boot_method == "auto") {
    method_code <- if (n <= 200L) 0L else if (n <= 5000L) 1L else 2L
  } else {
    method_code <- switch(boot_method,
      bca = 0L, percentile = 1L, parametric = 2L
    )
  }

  raw <- cpp_scale_ensemble_ci(x, n_boot, level, method_code)

  # Analytical CIs for each component estimator
  estimator_names <- c("sd_c4", "gmd", "mad_scaled", "iqr_scaled",
                       "sn", "qn", "robScale")

  analytical_lower <- numeric(7L)
  analytical_upper <- numeric(7L)
  for (j in seq_len(7L)) {
    if (j == 1L) {
      ci_obj <- .chisq_ci(raw$estimates[j], n, level)
    } else {
      ci_obj <- .analytical_ci(raw$estimates[j], n,
                               .are_values[[estimator_names[j]]],
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
