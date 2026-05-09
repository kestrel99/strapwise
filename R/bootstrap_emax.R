#' Bootstrap a fitted logistic Emax model
#'
#' Runs case-resampling bootstrap on a fitted `logistic_emax` object, returning
#' a matrix of bootstrapped parameter vectors alongside the original fit. The
#' result can be inspected via [print.boot_emax()] and [summary.boot_emax()].
#'
#' Each replicate resamples the rows of the original data (including all
#' covariate matrices) with replacement and refits [fit_logistic_emax()] using
#' the original estimates as warm-start values. Replicates that error, produce
#' warnings, or fail to converge are skipped and counted. A message is emitted
#' if any are skipped.
#'
#' @param fit A fitted [fit_logistic_emax()] object.
#' @param n_boot Number of bootstrap replicates. Default `1000L`.
#' @param seed Integer random seed passed to [set.seed()]. Default `42L`.
#' @param conf_level Confidence level for percentile intervals, strictly
#'   between 0 and 1. Default `0.95`.
#'
#' @return An S3 object of class `"boot_emax"`, a named list with:
#'   \describe{
#'     \item{`fit`}{The original `logistic_emax` object.}
#'     \item{`boot_theta`}{Numeric matrix of bootstrapped parameter vectors
#'       (`n_success` x `n_par`) in the optimisation parameterisation
#'       `(e0, emax, log_ec50, ...)`  -- the same scale as `fit$theta_hat`.}
#'     \item{`conf_level`}{The requested confidence level.}
#'     \item{`n_boot`}{Number of successful bootstrap replicates retained.}
#'   }
#'
#' @seealso [fit_logistic_emax()], [summary.boot_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#'
#' fit <- fit_logistic_emax(y, dose)
#' be  <- bootstrap_emax(fit, n_boot = 200, seed = 1)
#' print(be)
#' summary(be)
#'
#' @export
bootstrap_emax <- function(
  fit,
  n_boot     = 1000L,
  seed       = 42L,
  conf_level = 0.95
) {
  if (!inherits(fit, "logistic_emax"))
    stop("`fit` must be a fitted logistic_emax object.")
  if (conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be strictly between 0 and 1.")

  n         <- length(fit$x)
  n_par     <- fit$n_par
  th_names  <- names(fit$theta_hat)

  # Warm-start vector: convert log_ec50 back to natural scale so
  # fit_logistic_emax's start-handling can round-trip it correctly.
  warm <- as.numeric(fit$theta_hat)
  warm[fit$idx$log_ec50] <- exp(warm[fit$idx$log_ec50])

  # Pre-extract covariate matrices; NULL when 0-column (no covariates).
  Z_lin  <- if (ncol(fit$Z_lin)  > 0) fit$Z_lin  else NULL
  Z_emax <- if (ncol(fit$Z_emax) > 0) fit$Z_emax else NULL
  Z_ec50 <- if (ncol(fit$Z_ec50) > 0) fit$Z_ec50 else NULL

  boot_theta <- matrix(
    NA_real_,
    nrow     = n_boot,
    ncol     = n_par,
    dimnames = list(NULL, th_names)
  )
  n_failed <- 0L

  set.seed(seed)
  for (i in seq_len(n_boot)) {
    rows <- sample.int(n, n, replace = TRUE)

    warned   <- FALSE
    boot_fit <- tryCatch(
      withCallingHandlers(
        fit_logistic_emax(
          y           = fit$y[rows],
          x           = fit$x[rows],
          linear_covs = if (!is.null(Z_lin))  Z_lin[rows, , drop = FALSE]  else NULL,
          emax_covs   = if (!is.null(Z_emax)) Z_emax[rows, , drop = FALSE] else NULL,
          ec50_covs   = if (!is.null(Z_ec50)) Z_ec50[rows, , drop = FALSE] else NULL,
          start       = warm,
          hessian     = FALSE,
          control     = list(maxit = 500)
        ),
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

    boot_theta[i, ] <- boot_fit$theta_hat
  }

  if (n_failed > 0L)
    message(sprintf(
      "bootstrap_emax: %d replicate(s) skipped due to errors, warnings, or non-convergence.",
      n_failed
    ))

  boot_theta <- boot_theta[stats::complete.cases(boot_theta), , drop = FALSE]
  n_success  <- nrow(boot_theta)

  if (n_success < 100L)
    warning(sprintf(
      "Only %d bootstrap replicates succeeded; CIs may be unreliable.", n_success
    ))

  structure(
    list(
      fit        = fit,
      boot_theta = boot_theta,
      conf_level = conf_level,
      n_boot     = n_success
    ),
    class = "boot_emax"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `boot_emax` object
#'
#' @param x A `boot_emax` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [bootstrap_emax()], [summary.boot_emax()]
#'
#' @export
print.boot_emax <- function(x, ...) {
  cat("-- Bootstrap Logistic Emax -----------------------------------------\n")
  cat(sprintf("  N obs      : %d\n",  length(x$fit$x)))
  cat(sprintf("  Parameters : %d\n",  x$fit$n_par))
  cat(sprintf("  Replicates : %d successful\n", x$n_boot))
  cat(sprintf("  CI level   : %.0f%%\n", x$conf_level * 100))
  cat("\nOriginal-data estimates:\n")
  print(x$fit$coefficients)
  invisible(x)
}


#' Summarise a `boot_emax` object
#'
#' Returns and prints a data frame of original-data estimates alongside
#' bootstrapped standard errors, bias, and percentile confidence intervals.
#' The ec50 column of `boot_theta` is back-transformed from log scale before
#' computing summaries so all CIs are on the reporting scale.
#'
#' @param object A `boot_emax` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return A data frame (invisibly) with columns `Estimate`, `Boot_SE`,
#'   `Boot_Bias`, `CI_<level>%_lo`, and `CI_<level>%_hi`.
#'
#' @seealso [bootstrap_emax()], [print.boot_emax()]
#'
#' @export
summary.boot_emax <- function(object, ...) {
  alpha  <- 1 - object$conf_level
  ci_pct <- round(object$conf_level * 100)
  fit    <- object$fit

  # Convert boot_theta to reporting scale: back-transform the ec50 column.
  boot_rep <- object$boot_theta
  boot_rep[, fit$idx$log_ec50] <- exp(boot_rep[, fit$idx$log_ec50])
  colnames(boot_rep) <- names(fit$coefficients)

  coef_orig  <- fit$coefficients
  ci_lo      <- apply(boot_rep, 2, stats::quantile, probs = alpha / 2,       na.rm = TRUE)
  ci_hi      <- apply(boot_rep, 2, stats::quantile, probs = 1 - alpha / 2,   na.rm = TRUE)
  boot_se    <- apply(boot_rep, 2, stats::sd, na.rm = TRUE)
  boot_bias  <- colMeans(boot_rep, na.rm = TRUE) - coef_orig

  tbl <- data.frame(
    Estimate  = round(coef_orig,  4),
    Boot_SE   = round(boot_se,    4),
    Boot_Bias = round(boot_bias,  4),
    CI_lo     = round(ci_lo,      4),
    CI_hi     = round(ci_hi,      4),
    check.names = FALSE
  )
  names(tbl)[4:5] <- c(
    sprintf("CI_%d%%_lo", ci_pct),
    sprintf("CI_%d%%_hi", ci_pct)
  )

  cat("-- Bootstrap Emax Summary ------------------------------------------\n")
  cat(sprintf(
    "  Replicates : %d  |  CI level: %.0f%%\n\n",
    object$n_boot, object$conf_level * 100
  ))
  print(tbl)
  invisible(tbl)
}
