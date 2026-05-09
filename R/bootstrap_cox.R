#' Bootstrap a fitted Cox proportional hazards model
#'
#' Runs case-resampling bootstrap on a fitted [survival::coxph()] object,
#' returning a matrix of bootstrapped log hazard-ratio estimates alongside the
#' original fit. The result can be inspected via [print.boot_cox()] and
#' [summary.boot_cox()].
#'
#' @details
#' The function recovers the original data from `fit$call$data` so that the
#' full (untransformed) data frame is available for row resampling. If the call
#' environment cannot be evaluated it falls back to `model.frame(fit)` with a
#' warning. The `ties` method of the original fit is preserved in every
#' bootstrap refit.
#'
#' Replicates that error or produce warnings during refitting (e.g. Firth
#' correction triggers, collinearity warnings) are skipped and counted. If any
#' are skipped a message is emitted. Replicates with any `NA` coefficient are
#' dropped before the result is returned.
#'
#' @param fit A fitted [survival::coxph()] object.
#' @param n_boot Number of bootstrap replicates. Default `1000L`.
#' @param seed Integer random seed passed to [set.seed()]. Default `42L`.
#' @param conf_level Confidence level for percentile intervals, strictly
#'   between 0 and 1. Default `0.95`.
#'
#' @return An S3 object of class `"boot_cox"`, a named list with:
#'   \describe{
#'     \item{`fit`}{The original `coxph` object.}
#'     \item{`boot_coefs`}{Numeric matrix of bootstrapped log hazard-ratio
#'       estimates (`n_success` x `p`), after dropping incomplete rows.}
#'     \item{`data`}{The original data frame recovered from the fit call.}
#'     \item{`formula`}{The model formula.}
#'     \item{`conf_level`}{The requested confidence level.}
#'     \item{`n_boot`}{Number of successful bootstrap replicates retained.}
#'   }
#'
#' @seealso [summary.boot_cox()], [print.boot_cox()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L   # recode to 0/1
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex + ph.ecog,
#'     data = lung
#'   )
#'   bc <- bootstrap_cox(fit, n_boot = 200, seed = 1)
#'   print(bc)
#'   summary(bc)
#' }
#'
#' @export
bootstrap_cox <- function(
  fit,
  n_boot     = 1000L,
  seed       = 42L,
  conf_level = 0.95
) {
  if (!inherits(fit, "coxph"))
    stop("`fit` must be a fitted coxph object.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1.")

  formula    <- formula(fit)
  ties_method <- fit$method   # "efron", "breslow", or "exact"

  # ---- Recover original data -------------------------------------------------
  data <- tryCatch(
    eval(fit$call$data, envir = environment(formula)),
    error = function(e) {
      warning(
        "bootstrap_cox: could not recover original data from fit$call$data; ",
        "falling back to model.frame(fit)."
      )
      stats::model.frame(fit)
    }
  )

  n          <- nrow(data)
  coef_names <- names(coef(fit))

  # ---- Bootstrap loop --------------------------------------------------------
  set.seed(seed)
  boot_coefs <- matrix(
    NA_real_,
    nrow     = n_boot,
    ncol     = length(coef_names),
    dimnames = list(NULL, coef_names)
  )
  n_failed <- 0L

  for (i in seq_len(n_boot)) {
    rows      <- sample.int(n, n, replace = TRUE)
    boot_data <- data[rows, , drop = FALSE]

    warned   <- FALSE
    boot_fit <- tryCatch(
      withCallingHandlers(
        survival::coxph(formula, data = boot_data, ties = ties_method),
        warning = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )

    if (is.null(boot_fit) || warned) {
      n_failed <- n_failed + 1L
      next
    }

    # Name-match guards against dropped terms in degenerate resamples
    matched <- intersect(names(coef(boot_fit)), coef_names)
    boot_coefs[i, matched] <- coef(boot_fit)[matched]
  }

  if (n_failed > 0L)
    message(sprintf(
      "bootstrap_cox: %d replicate(s) skipped due to errors or warnings.",
      n_failed
    ))

  boot_coefs <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]
  n_success  <- nrow(boot_coefs)

  if (n_success < 100L)
    warning(sprintf(
      "Only %d bootstrap replicates succeeded; CIs may be unreliable.",
      n_success
    ))

  structure(
    list(
      fit        = fit,
      boot_coefs = boot_coefs,
      data       = data,
      formula    = formula,
      conf_level = conf_level,
      n_boot     = n_success
    ),
    class = "boot_cox"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `boot_cox` object
#'
#' Displays a compact summary of the bootstrap run and the original-data log
#' hazard-ratio estimates.
#'
#' @param x A `boot_cox` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [bootstrap_cox()], [summary.boot_cox()]
#'
#' @export
print.boot_cox <- function(x, ...) {
  cat("-- Bootstrap Cox Proportional Hazards ------------------------------\n")
  cat(sprintf("  Formula   : %s\n", deparse(x$formula)))
  cat(sprintf("  N obs     : %d\n", x$fit$n))
  cat(sprintf("  N events  : %d\n", x$fit$nevent))
  cat(sprintf("  Ties      : %s\n", x$fit$method))
  cat(sprintf("  Replicates: %d successful\n", x$n_boot))
  cat(sprintf("  CI level  : %.0f%%\n", x$conf_level * 100))
  cat("\nOriginal-data log hazard ratios:\n")
  print(coef(x$fit))
  invisible(x)
}


#' Summarise a `boot_cox` object
#'
#' Returns and prints a data frame of original-data estimates alongside
#' bootstrapped standard errors, bias, and percentile confidence intervals.
#' Results are presented on both the log hazard-ratio scale (for inference)
#' and the hazard-ratio scale (for interpretation).
#'
#' @param object A `boot_cox` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return A data frame (invisibly) with columns `log_HR`, `HR`, `Boot_SE`,
#'   `Boot_Bias`, `CI_<level>%_lo`, `CI_<level>%_hi`, `HR_<level>%_lo`, and
#'   `HR_<level>%_hi`, where `<level>` is `round(object$conf_level * 100)`.
#'
#' @seealso [bootstrap_cox()], [print.boot_cox()]
#'
#' @export
summary.boot_cox <- function(object, ...) {
  alpha      <- 1 - object$conf_level
  ci_pct     <- round(object$conf_level * 100)
  coef_orig  <- coef(object$fit)

  ci_lo    <- apply(object$boot_coefs, 2, stats::quantile,
                    probs = alpha / 2,       na.rm = TRUE)
  ci_hi    <- apply(object$boot_coefs, 2, stats::quantile,
                    probs = 1 - alpha / 2,   na.rm = TRUE)
  boot_se   <- apply(object$boot_coefs, 2, stats::sd,   na.rm = TRUE)
  boot_bias <- colMeans(object$boot_coefs, na.rm = TRUE) - coef_orig

  tbl <- data.frame(
    log_HR    = round(coef_orig,        4),
    HR        = round(exp(coef_orig),   4),
    Boot_SE   = round(boot_se,          4),
    Boot_Bias = round(boot_bias,        4),
    CI_lo     = round(ci_lo,            4),
    CI_hi     = round(ci_hi,            4),
    HR_lo     = round(exp(ci_lo),       4),
    HR_hi     = round(exp(ci_hi),       4),
    check.names = FALSE
  )
  names(tbl)[5:8] <- c(
    sprintf("CI_%d%%_lo",    ci_pct),
    sprintf("CI_%d%%_hi",    ci_pct),
    sprintf("HR_%d%%_lo",    ci_pct),
    sprintf("HR_%d%%_hi",    ci_pct)
  )

  cat("-- Bootstrap Cox Summary -------------------------------------------\n")
  cat(sprintf(
    "  Replicates : %d  |  CI level: %.0f%%\n\n",
    object$n_boot, object$conf_level * 100
  ))
  print(tbl)
  invisible(tbl)
}
