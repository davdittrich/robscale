# Not exported; called by individual estimator functions when ci = TRUE.
# Source: Rousseeuw & Croux (1993) and standard references.
.are_values <- c(
  sd_c4      = 1.00,
  gmd        = 0.98,
  mad_scaled = 0.368,
  iqr_scaled = 0.37,
  sn         = 0.58,
  qn         = 0.82,
  robScale   = 0.55,
  adm        = 0.88
)

# ARE-based normal approximation CI (all estimators except sd_c4)
.analytical_ci <- function(estimate, n, are, level, method) {
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  se <- estimate / sqrt(2 * n * are)
  structure(
    list(
      estimate    = estimate,
      ci          = c(lower = estimate - z * se, upper = estimate + z * se),
      level       = level,
      method      = method,
      boot_method = "analytical"
    ),
    class = "robscale_ci"
  )
}

# Exact chi-squared CI for sd_c4
.chisq_ci <- function(estimate, n, level) {
  alpha <- 1 - level
  df <- n - 1L
  structure(
    list(
      estimate    = estimate,
      ci          = c(
        lower = estimate * sqrt(df / qchisq(1 - alpha / 2, df)),
        upper = estimate * sqrt(df / qchisq(alpha / 2, df))
      ),
      level       = level,
      method      = "sd_c4",
      boot_method = "analytical"
    ),
    class = "robscale_ci"
  )
}

#' Print a robscale_ci object
#'
#' @param x A \code{"robscale_ci"} object.
#' @param digits Number of significant digits to display. Default: 4.
#' @param ... Further arguments passed to or from other methods.
#' @return The object \code{x}, invisibly.
#' @export
print.robscale_ci <- function(x, digits = 4, ...) {
  cat(sprintf("%s estimate: %s\n", x$method,
              formatC(x$estimate, digits = digits, format = "f")))
  bm <- if (!is.null(x$boot_method)) x$boot_method else "analytical"
  cat(sprintf("%g%% CI (%s): [%s, %s]\n",
    x$level * 100,
    bm,
    formatC(x$ci[["lower"]], digits = digits, format = "f"),
    formatC(x$ci[["upper"]], digits = digits, format = "f")))
  invisible(x)
}

#' Print a robscale_ensemble_ci object
#'
#' @param x A \code{"robscale_ensemble_ci"} object.
#' @param digits Number of significant digits to display. Default: 4.
#' @param ... Further arguments passed to or from other methods.
#' @return The object \code{x}, invisibly.
#' @export
print.robscale_ensemble_ci <- function(x, digits = 4, ...) {
  cat(sprintf("Ensemble estimate: %s\n",
              formatC(x$estimate, digits = digits, format = "f")))
  cat(sprintf("%g%% bootstrap CI (%s): [%s, %s]\n\n",
    x$level * 100,
    x$boot_method,
    formatC(x$ci[["lower"]], digits = digits, format = "f"),
    formatC(x$ci[["upper"]], digits = digits, format = "f")))
  cat("Component estimators:\n")
  print(x$estimators, digits = digits, row.names = FALSE)
  invisible(x)
}
