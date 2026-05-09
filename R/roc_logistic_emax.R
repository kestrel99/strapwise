#' ROC curve for a logistic Emax model
#'
#' Computes sensitivity and specificity across a threshold grid, plots the ROC
#' curve with an optional shaded AUC area, and marks the Youden-optimal
#' threshold. Returns an object of class `"roc_logistic_emax"` that holds the
#' plot, the ROC data frame, the AUC, and the best-threshold row.
#'
#' AUC is computed via the trapezoidal rule. The confidence interval (when
#' `ci = TRUE`) uses the Hanley-McNeil approximation.
#'
#' @param object A `logistic_emax` object returned by [fit_logistic_emax()].
#' @param thresholds Numeric vector of classification thresholds. When `NULL`
#'   (default), 501 evenly spaced values from 1 to 0 are used.
#' @param ci Logical; if `TRUE` (default) a 95% Hanley-McNeil confidence
#'   interval is appended to the AUC subtitle.
#' @param ref_line Logical; if `TRUE` (default) a dashed diagonal chance line
#'   is drawn.
#' @param color Line colour. Default `"steelblue"`.
#' @param fill Fill colour for the AUC area. Default `"steelblue"`.
#' @param alpha Opacity of the AUC ribbon. Default `0.2`.
#' @param title Character string for the plot title. Default
#'   `"ROC Curve - Logistic Emax Model"`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return An S3 object of class `"roc_logistic_emax"`, a named list with:
#'   \describe{
#'     \item{`plot`}{A `ggplot` object.}
#'     \item{`roc_df`}{Data frame with columns `threshold`, `sensitivity`,
#'       `specificity`, `fpr`, and `youden`.}
#'     \item{`auc`}{Scalar AUC.}
#'     \item{`best_threshold`}{Single-row data frame for the Youden-optimal
#'       threshold.}
#'   }
#'
#' @seealso [fit_logistic_emax()], [print.roc_logistic_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#' fit  <- fit_logistic_emax(y, dose)
#'
#' roc <- roc_logistic_emax(fit)
#' print(roc)
#'
#' @export
roc_logistic_emax <- function(
  object,
  thresholds = NULL,
  ci         = TRUE,
  ref_line   = TRUE,
  color      = "steelblue",
  fill       = "steelblue",
  alpha      = 0.2,
  title      = "ROC Curve - Logistic Emax Model",
  ...
) {
  y_obs  <- object$y
  y_pred <- object$fitted.values

  if (is.null(thresholds)) {
    thresholds <- seq(1, 0, length.out = 501)
  }

  roc_df <- do.call(rbind, lapply(thresholds, function(thr) {
    pred_pos <- y_pred >= thr
    tp <- sum( pred_pos &  y_obs == 1)
    fp <- sum( pred_pos &  y_obs == 0)
    fn <- sum(!pred_pos &  y_obs == 1)
    tn <- sum(!pred_pos &  y_obs == 0)
    data.frame(
      threshold   = thr,
      sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
      specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    )
  }))
  roc_df$fpr    <- 1 - roc_df$specificity
  roc_df        <- roc_df[order(roc_df$fpr, roc_df$sensitivity), ]

  auc <- with(roc_df, {
    idx <- order(fpr)
    x   <- fpr[idx]
    y   <- sensitivity[idx]
    sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
  })

  if (ci) {
    n1     <- sum(y_obs == 1)
    n0     <- sum(y_obs == 0)
    q1     <- auc / (2 - auc)
    q2     <- 2 * auc^2 / (1 + auc)
    se_auc <- sqrt(
      (auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) /
      (n1 * n0)
    )
    z <- stats::qnorm(0.975)
    auc_str <- sprintf(
      "\nAUC = %.3f (95%% CI: %.3f-%.3f)",
      auc, pmax(0, auc - z * se_auc), pmin(1, auc + z * se_auc)
    )
  } else {
    auc_str <- sprintf("\nAUC = %.3f", auc)
  }

  roc_df$youden <- roc_df$sensitivity + roc_df$specificity - 1
  best <- roc_df[which.max(roc_df$youden), ]

  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = sensitivity)) +
    ggplot2::geom_area(fill = fill, alpha = alpha) +
    ggplot2::geom_line(color = color, linewidth = 1.1) +
    ggplot2::geom_point(
      data  = best,
      ggplot2::aes(x = fpr, y = sensitivity),
      color = "firebrick", size = 3, shape = 16
    ) +
    ggplot2::geom_label(
      data = best,
      ggplot2::aes(
        x     = fpr + 0.03,
        y     = sensitivity - 0.05,
        label = sprintf(
          "Thr = %.2f\nSe = %.2f, Sp = %.2f",
          threshold, sensitivity, specificity
        )
      ),
      hjust = 0, size = 3.2, color = "firebrick",
      fill = "white", label.size = 0.3
    ) +
    { if (ref_line)
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             linetype = "dashed", color = "grey50")
    } +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(
      x        = "1 - Specificity (FPR)",
      y        = "Sensitivity (TPR)",
      title    = title,
      subtitle = paste0(
        sprintf(
          "E0 = %.2f | Emax = %.2f | EC50 = %.2f",
          object$coefficients[1],
          object$coefficients[2],
          object$coefficients[3]
        ),
        auc_str
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(aspect.ratio = 1)

  invisible(structure(
    list(plot = p, roc_df = roc_df, auc = auc, best_threshold = best),
    class = "roc_logistic_emax"
  ))
}


#' Print a `roc_logistic_emax` object
#'
#' Prints the AUC and optimal threshold, then renders the ROC plot.
#'
#' @param x A `roc_logistic_emax` object.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [roc_logistic_emax()]
#'
#' @export
print.roc_logistic_emax <- function(x, ...) {
  cat(sprintf("AUC: %.4f\n", x$auc))
  cat(sprintf(
    "Optimal threshold (Youden): %.3f  Se = %.3f  Sp = %.3f\n",
    x$best_threshold$threshold,
    x$best_threshold$sensitivity,
    x$best_threshold$specificity
  ))
  print(x$plot)
  invisible(x)
}
