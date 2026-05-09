#' Backwards stepwise model reduction by likelihood ratio test
#'
#' Iteratively removes the least significant term from a model until all
#' remaining terms satisfy \code{p < alpha}, or no testable terms remain.
#' Supports binomial [glm()] and [survival::coxph()] objects via a unified
#' interface.
#'
#' At each step every candidate term is tested by comparing the current model
#' to the model with that term omitted (LRT):
#' \deqn{\chi^2 = -2\bigl[\ell(\text{reduced}) - \ell(\text{current})\bigr]}
#' with degrees of freedom equal to the difference in the number of estimated
#' coefficients.  Multi-level factors are handled correctly because the
#' coefficient count correctly captures their effective degrees of freedom.
#' The term with the largest p-value \eqn{\geq \alpha} is removed.  Iteration
#' stops when all remaining (testable) terms satisfy \eqn{p < \alpha}.
#'
#' The function recovers the original data from \code{fit\$call\$data} so that
#' each reduced model is refit on identical data.  If recovery fails it falls
#' back to [model.frame()] with a message.
#'
#' @param fit A fitted [glm()] with `family = binomial` (any link), or a
#'   fitted [survival::coxph()] object.
#' @param alpha Significance threshold.  Terms with LRT p-value
#'   \eqn{\geq \alpha} are candidates for removal.  Default `0.05`.
#' @param keep_vars Optional character vector of term names that must not be
#'   dropped regardless of their p-value (e.g. treatment, strata variables).
#' @param verbose Logical; print a step-by-step trace. Default `TRUE`.
#' @param max_iterations Maximum number of steps before stopping with a
#'   warning. Default `50L`.
#'
#' @return An S3 object of class `"stepwise_reduction"`, a named list with:
#'   \describe{
#'     \item{`final_model`}{The model after reduction.}
#'     \item{`original_model`}{The model passed as `fit`.}
#'     \item{`step_table`}{Data frame with one row per dropped term and
#'       columns: `step`, `term`, `df`, `AIC`, `dAIC`, `neg2LL`,
#'       `d_neg2LL` (LRT statistic), `p`, `n_terms` (terms remaining).}
#'     \item{`dropped`}{Character vector of removed terms, in drop order.}
#'     \item{`kept`}{Value of `keep_vars`.}
#'     \item{`alpha`}{Significance threshold used.}
#'     \item{`model_type`}{`"Logistic"` or `"Cox"`.}
#'   }
#'
#' @seealso [print.stepwise_reduction()], [summary.stepwise_reduction()],
#'   [bootstrap_logistic()], [bootstrap_cox()]
#'
#' @examples
#' # ---- Logistic regression ---------------------------------------------------
#' set.seed(1)
#' df <- data.frame(
#'   y   = rbinom(200, 1, 0.35),
#'   age = rnorm(200, 50, 10),
#'   bmi = rnorm(200, 27,  5),
#'   sex = factor(sample(c("F", "M"), 200, replace = TRUE))
#' )
#' fit_log <- glm(y ~ age + bmi + sex, data = df, family = binomial)
#' sr      <- stepwise_reduce(fit_log, alpha = 0.10)
#' print(sr)
#'
#' # ---- Cox proportional hazards ----------------------------------------------
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   lung_cc <- lung[complete.cases(
#'     lung[, c("time", "status", "age", "sex", "ph.ecog", "wt.loss")]
#'   ), ]
#'   fit_cox <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex + ph.ecog + wt.loss,
#'     data = lung_cc
#'   )
#'   sr_cox <- stepwise_reduce(fit_cox, alpha = 0.05)
#'   summary(sr_cox)
#' }
#'
#' @export
stepwise_reduce <- function(
  fit,
  alpha          = 0.05,
  keep_vars      = NULL,
  verbose        = TRUE,
  max_iterations = 50L
) {
  # ---- Validate ----------------------------------------------------------------
  is_cox <- inherits(fit, "coxph")
  is_glm <- inherits(fit, "glm") &&
              identical(fit$family$family, "binomial")

  if (!is_cox && !is_glm)
    stop("`fit` must be a binomial glm or a coxph object.")
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be strictly between 0 and 1.")
  if (!is.null(keep_vars)) {
    if (!is.character(keep_vars))
      stop("`keep_vars` must be a character vector or NULL.")
    unknown <- setdiff(keep_vars, attr(stats::terms(fit), "term.labels"))
    if (length(unknown) > 0L)
      warning(sprintf(
        "keep_vars contains term(s) not in the model: %s",
        paste(unknown, collapse = ", ")
      ))
  }

  model_type  <- if (is_cox) "Cox" else "Logistic"
  ties_method <- if (is_cox) fit$method else NULL

  # ---- Recover training data ---------------------------------------------------
  fit_data <- tryCatch(
    eval(fit$call$data, envir = environment(stats::formula(fit))),
    error = function(e) {
      message("stepwise_reduce: cannot recover data from fit$call$data; ",
              "falling back to model.frame(fit).")
      stats::model.frame(fit)
    }
  )

  # ---- Hierarchy filter: terms that can be dropped without violating order ----
  # A term is droppable only when none of the remaining terms is a higher-order
  # interaction that contains every predictor of that term.  This mirrors the
  # behaviour of drop1() (which never tests a main effect while its interaction
  # is still in the model).
  .droppable_terms <- function(m) {
    trm    <- stats::terms(m)
    labels <- attr(trm, "term.labels")
    if (length(labels) == 0L) return(character(0L))
    fac <- attr(trm, "factors")   # predictors x terms matrix (0/1)
    ord <- attr(trm, "order")     # interaction order per term
    vapply(seq_along(labels), function(i) {
      hi <- which(ord > ord[i])
      if (length(hi) == 0L) return(TRUE)
      pred_i <- which(fac[, i] > 0)
      # blocked if any higher-order term j contains ALL predictors of term i
      !any(vapply(hi,
                  function(j) all(fac[pred_i, j] > 0),
                  logical(1L)))
    }, logical(1L)) |> (\(keep) labels[keep])()
  }

  # ---- LRT helper: test removing one term from `current` ----------------------
  # Captures fit_data, is_cox, ties_method from enclosing scope.
  .lrt_one <- function(current, term) {
    new_f <- stats::update.formula(
      stats::formula(current),
      stats::as.formula(paste(". ~ . -", term))
    )

    reduced <- if (is_cox) {
      tryCatch(
        withCallingHandlers(
          survival::coxph(new_f, data = fit_data, ties = ties_method),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) NULL
      )
    } else {
      tryCatch(
        withCallingHandlers(
          stats::glm(new_f, data = fit_data, family = current$family),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) NULL
      )
    }

    if (is.null(reduced))
      return(list(p = NA_real_, df = NA_integer_, lrt = NA_real_,
                  model = NULL))

    df_diff  <- max(1L,
                    length(stats::coef(current)) -
                    length(stats::coef(reduced)))
    lrt_stat <- max(0,
                    -2 * as.numeric(stats::logLik(reduced) -
                                    stats::logLik(current)))
    p_val    <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)

    list(p = p_val, df = df_diff, lrt = lrt_stat, model = reduced)
  }

  # ---- Main loop ---------------------------------------------------------------
  current    <- fit
  dropped    <- character(0L)
  step_rows  <- list()
  converged  <- FALSE

  if (verbose) {
    pad <- strrep("-", max(0, 48 - nchar(model_type)))
    cat(sprintf("-- Stepwise %s model reduction  (alpha = %g) %s\n",
                model_type, alpha, pad))
    cat(sprintf("   Starting terms : %d\n",
                length(attr(stats::terms(fit), "term.labels"))))
    if (!is.null(keep_vars) && length(keep_vars) > 0L)
      cat(sprintf("   Protected      : %s\n",
                  paste(keep_vars, collapse = ", ")))
    cat(strrep("-", 60), "\n")
  }

  for (iter in seq_len(max_iterations)) {
    term_labels <- attr(stats::terms(current), "term.labels")
    testable    <- setdiff(.droppable_terms(current), keep_vars)

    if (length(testable) == 0L) {
      if (verbose) cat("   No testable terms remain. Stopping.\n")
      converged <- TRUE
      break
    }

    prev_aic    <- AIC(current)
    prev_neg2ll <- -2 * as.numeric(stats::logLik(current))

    # LRT for all candidate terms
    tests  <- lapply(stats::setNames(testable, testable),
                     function(trm) .lrt_one(current, trm))
    p_vals <- vapply(tests, `[[`, numeric(1L), "p")

    # Failed refits (NA p) are treated as significant  -- keep them
    p_vals[is.na(p_vals)] <- 0

    droppable <- p_vals[p_vals >= alpha]

    if (length(droppable) == 0L) {
      if (verbose)
        cat(sprintf(
          "   All remaining terms p < %g. Final model has %d term(s).\n",
          alpha, length(term_labels)
        ))
      converged <- TRUE
      break
    }

    drop_term <- names(which.max(droppable))
    drop_res  <- tests[[drop_term]]

    current <- drop_res$model
    dropped <- c(dropped, drop_term)

    cur_aic    <- AIC(current)
    cur_neg2ll <- -2 * as.numeric(stats::logLik(current))
    n_left     <- length(attr(stats::terms(current), "term.labels"))

    step_rows[[iter]] <- data.frame(
      step     = iter,
      term     = drop_term,
      df       = drop_res$df,
      AIC      = cur_aic,
      dAIC     = cur_aic    - prev_aic,
      neg2LL   = cur_neg2ll,
      d_neg2LL = drop_res$lrt,
      p        = drop_res$p,
      n_terms  = n_left,
      stringsAsFactors = FALSE
    )

    if (verbose)
      cat(sprintf(
        "   Step %2d  %-25s  df=%d  LRT=%7.3f  p=%s\n",
        iter, drop_term, drop_res$df, drop_res$lrt,
        format.pval(drop_res$p, digits = 3, eps = 0.001)
      ))

    if (n_left == 0L) {
      converged <- TRUE
      break
    }
  }

  if (!converged)
    warning(sprintf(
      "stepwise_reduce: max_iterations = %d reached without convergence.",
      max_iterations
    ))

  step_table <- if (length(step_rows) > 0L) {
    out <- do.call(rbind, step_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      step = integer(0), term = character(0), df = integer(0),
      AIC = numeric(0), dAIC = numeric(0), neg2LL = numeric(0),
      d_neg2LL = numeric(0), p = numeric(0), n_terms = integer(0),
      stringsAsFactors = FALSE
    )
  }

  if (verbose) {
    cat(strrep("-", 60), "\n")
    n_final <- length(attr(stats::terms(current), "term.labels"))
    cat(sprintf("   Dropped %d term(s); %d term(s) remain.\n",
                length(dropped), n_final))
  }

  structure(
    list(
      final_model    = current,
      original_model = fit,
      step_table     = step_table,
      dropped        = dropped,
      kept           = keep_vars,
      alpha          = alpha,
      model_type     = model_type
    ),
    class = "stepwise_reduction"
  )
}


# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' Print a `stepwise_reduction` object
#'
#' Displays the reduction trace (step table) in a compact, human-readable
#' format.
#'
#' @param x A `stepwise_reduction` object returned by [stepwise_reduce()].
#' @param digits Integer; number of decimal places for numeric columns.
#'   Default `4`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [stepwise_reduce()], [summary.stepwise_reduction()]
#'
#' @export
print.stepwise_reduction <- function(x, digits = 4, ...) {
  n_orig  <- length(attr(stats::terms(x$original_model), "term.labels"))
  n_final <- length(attr(stats::terms(x$final_model),    "term.labels"))
  n_steps <- if (is.null(x$step_table) || nrow(x$step_table) == 0L) 0L
             else nrow(x$step_table)

  cat(sprintf("-- Stepwise %s Model Reduction\n", x$model_type))
  cat(sprintf("   alpha = %g  |  %d step(s)  |  terms: %d -> %d\n",
              x$alpha, n_steps, n_orig, n_final))

  if (!is.null(x$kept) && length(x$kept) > 0L)
    cat(sprintf("   Protected: %s\n", paste(x$kept, collapse = ", ")))

  if (n_steps == 0L) {
    cat("   No terms were dropped.\n")
    return(invisible(x))
  }

  cat(sprintf("   Dropped  : %s\n\n", paste(x$dropped, collapse = ", ")))

  tbl <- x$step_table
  tbl$AIC      <- round(tbl$AIC,      digits)
  tbl$dAIC     <- round(tbl$dAIC,     digits)
  tbl$neg2LL   <- round(tbl$neg2LL,   digits)
  tbl$d_neg2LL <- round(tbl$d_neg2LL, digits)
  tbl$p        <- format.pval(tbl$p, digits = max(2L, digits - 1L),
                               eps = 0.001)

  names(tbl)[names(tbl) == "neg2LL"]   <- "-2LL"
  names(tbl)[names(tbl) == "d_neg2LL"] <- "d(-2LL)"
  names(tbl)[names(tbl) == "n_terms"]  <- "terms_left"

  print(tbl, row.names = FALSE)
  invisible(x)
}


#' Summarise a `stepwise_reduction` object
#'
#' Prints the step table followed by a side-by-side comparison of the original
#' and final model fit statistics, and an overall LRT between the two models.
#'
#' @param object A `stepwise_reduction` object returned by [stepwise_reduce()].
#' @param digits Integer; number of decimal places for numeric columns.
#'   Default `4`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `object`, invisibly.
#'
#' @seealso [stepwise_reduce()], [print.stepwise_reduction()]
#'
#' @export
summary.stepwise_reduction <- function(object, digits = 4, ...) {
  print(object, digits = digits)

  orig  <- object$original_model
  final <- object$final_model

  cat("\n-- Model Comparison (Original vs Final) ", strrep("-", 20), "\n")
  cat(sprintf("   Original : %s\n", deparse(stats::formula(orig))))
  cat(sprintf("   Final    : %s\n", deparse(stats::formula(final))))

  # ---- Fit statistics ----------------------------------------------------------
  aic_o  <- AIC(orig);        aic_f  <- AIC(final)
  ll_o   <- -2 * as.numeric(stats::logLik(orig))
  ll_f   <- -2 * as.numeric(stats::logLik(final))
  bic_o  <- tryCatch(BIC(orig),  error = function(e) NA_real_)
  bic_f  <- tryCatch(BIC(final), error = function(e) NA_real_)

  cat(sprintf(
    "\n   %-12s  %10s  %10s  %10s\n",
    "", "Original", "Final", "Delta"
  ))
  cat(sprintf("   %-12s  %10.3f  %10.3f  %+10.3f\n",
              "AIC", aic_o, aic_f, aic_f - aic_o))
  if (!is.na(bic_o))
    cat(sprintf("   %-12s  %10.3f  %10.3f  %+10.3f\n",
                "BIC", bic_o, bic_f, bic_f - bic_o))
  cat(sprintf("   %-12s  %10.3f  %10.3f  %+10.3f\n",
              "-2 log-lik", ll_o, ll_f, ll_f - ll_o))

  # ---- Overall LRT (original vs final) ----------------------------------------
  df_overall  <- max(1L,
                     length(stats::coef(orig)) -
                     length(stats::coef(final)))
  lrt_overall <- ll_f - ll_o
  p_overall   <- stats::pchisq(lrt_overall, df = df_overall,
                                lower.tail = FALSE)
  cat(sprintf(
    "\n   Overall LRT (original vs final): chi2(%d) = %.3f, p = %s\n",
    df_overall, lrt_overall,
    format.pval(p_overall, digits = 3, eps = 0.001)
  ))

  # ---- Model-type-specific extras ----------------------------------------------
  if (object$model_type == "Cox") {
    cat(sprintf("   Events / N   : %d / %d\n", orig$nevent, orig$n))
    cat(sprintf("   Ties method  : %s\n", orig$method))
  } else {
    n_obs    <- nrow(stats::model.frame(orig))
    n_events <- sum(orig$y, na.rm = TRUE)
    cat(sprintf("   N obs / events: %d / %d\n", n_obs, n_events))
    cat(sprintf("   Residual df   : %d -> %d\n",
                stats::df.residual(orig), stats::df.residual(final)))
  }

  invisible(object)
}
