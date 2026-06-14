#' ROC curve(s) for one or more logistic Emax models
#'
#' Computes ROC curves for one or more `logistic_emax` objects, overlays them
#' on a single panel, and returns ROC data, AUC values, and the ggplot object.
#' AUC is computed via the trapezoidal rule; confidence intervals use the
#' Hanley-McNeil approximation.
#'
#' For a single model the plot matches the original behaviour: a shaded AUC
#' area, a subtitle with core parameters (E0, Emax, EC50) and AUC, and a
#' labelled Youden-optimal threshold point.  For multiple models the curves are
#' superimposed with a legend entry per model showing its AUC.
#'
#' @param object A `logistic_emax` object, **or** a named list of
#'   `logistic_emax` objects to compare.
#' @param labels Optional character vector of display labels, one per model.
#'   Defaults to `names(object)`; unnamed models are labelled `"Model 1"`,
#'   `"Model 2"`, etc.
#' @param thresholds Numeric vector of classification thresholds. When `NULL`
#'   (default), 501 evenly spaced values from 1 to 0 are used.
#' @param ci Logical; if `TRUE` (default) a 95% Hanley-McNeil confidence
#'   interval is appended to each AUC summary.
#' @param ref_line Logical; if `TRUE` (default) a dashed diagonal chance line
#'   is drawn.
#' @param fill Logical or `NULL`.  When `NULL` (default), area shading is
#'   applied for a single model and suppressed for multiple models.
#' @param alpha Opacity of the AUC shading. Default `0.2`.
#' @param colors Optional character vector of line colours, one per model.
#'   Defaults to the ggplot2 discrete colour scale.
#' @param linewidth Line width. Default `1.1`.
#' @param title Plot title. Default `"ROC Curve - Logistic Emax Model"`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return An S3 object of class `"roc_logistic_emax"`, a named list with:
#'   \describe{
#'     \item{`plot`}{A `ggplot` object.}
#'     \item{`roc_list`}{Named list of ROC data frames (one per model), each
#'       with columns `threshold`, `sensitivity`, `specificity`, `fpr`, and
#'       `youden`.}
#'     \item{`roc_df`}{Alias for `roc_list[[1]]`; retained for single-model
#'       backward compatibility.}
#'     \item{`auc`}{Named numeric vector of AUC values.}
#'     \item{`best_threshold`}{Named list of single-row data frames for the
#'       Youden-optimal threshold of each model.}
#'   }
#'
#' @seealso [fit_logistic_emax()], [print.roc_logistic_emax()],
#'   [roc_logistic()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#' fit  <- fit_logistic_emax(y, dose)
#'
#' # Single model
#' roc <- roc_logistic_emax(fit)
#' print(roc)
#'
#' # Two models superimposed
#' wt   <- rnorm(n, 70, 10)
#' fit2 <- fit_logistic_emax(y, dose, linear_covs = cbind(weight = wt))
#' roc2 <- roc_logistic_emax(list(Base = fit, WithWeight = fit2))
#' print(roc2)
#'
#' @export
roc_logistic_emax <- function(
  object,
  labels     = NULL,
  thresholds = NULL,
  ci         = TRUE,
  ref_line   = TRUE,
  fill       = NULL,
  alpha      = 0.2,
  colors     = NULL,
  linewidth  = 1.1,
  title      = "ROC Curve - Logistic Emax Model",
  ...
) {
  # ---- Normalise to a named list ---------------------------------------------
  if (inherits(object, "logistic_emax")) {
    models <- list(object)
  } else if (is.list(object)) {
    models <- object
  } else {
    stop("`object` must be a `logistic_emax` object or a named list of them.")
  }

  n_models <- length(models)
  if (n_models == 0L) stop("`object` must contain at least one model.")

  for (i in seq_len(n_models)) {
    if (!inherits(models[[i]], "logistic_emax"))
      stop("Element ", i, " of `object` is not a `logistic_emax` object.")
  }

  # ---- Labels ----------------------------------------------------------------
  if (is.null(labels)) labels <- names(models)
  if (is.null(labels)) labels <- character(n_models)
  unnamed <- is.na(labels) | nchar(labels) == 0L
  labels[unnamed] <- paste0("Model ", which(unnamed))
  if (length(labels) != n_models)
    stop("`labels` must have the same length as `object`.")
  if (anyDuplicated(labels))
    stop("`labels` must be unique.")

  # ---- Warn if outcomes differ across models ---------------------------------
  if (n_models > 1L) {
    y_ref <- models[[1L]]$y
    for (i in seq.int(2L, n_models)) {
      if (!identical(models[[i]]$y, y_ref))
        warning(
          "Model '", labels[i], "' has a different outcome vector than '",
          labels[1L], "'. AUC values may not be directly comparable."
        )
    }
  }

  if (is.null(thresholds))
    thresholds <- seq(1, 0, length.out = 501L)

  if (is.null(fill)) fill <- (n_models == 1L)

  if (!is.null(colors) && length(colors) < n_models)
    stop("`colors` must have at least ", n_models, " element(s).")

  # ---- ROC helpers -----------------------------------------------------------
  roc_for_model <- function(y_obs, y_pred) {
    df <- do.call(rbind, lapply(thresholds, function(thr) {
      pos <- y_pred >= thr
      tp  <- sum( pos &  y_obs == 1)
      fp  <- sum( pos &  y_obs == 0)
      fn  <- sum(!pos &  y_obs == 1)
      tn  <- sum(!pos &  y_obs == 0)
      data.frame(
        threshold   = thr,
        sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
        specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      )
    }))
    df$fpr    <- 1 - df$specificity
    df        <- df[order(df$fpr, df$sensitivity), ]
    df$youden <- df$sensitivity + df$specificity - 1
    df
  }

  trapz_auc <- function(df) {
    idx <- order(df$fpr)
    x   <- df$fpr[idx]
    y   <- df$sensitivity[idx]
    sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
  }

  hm_ci <- function(auc, n1, n0) {
    q1  <- auc / (2 - auc)
    q2  <- 2 * auc^2 / (1 + auc)
    num <- auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)
    se  <- sqrt(num / (n1 * n0))
    z   <- stats::qnorm(0.975)
    c(lo = pmax(0, auc - z * se), hi = pmin(1, auc + z * se))
  }

  # ---- Compute ROC per model -------------------------------------------------
  roc_list      <- setNames(vector("list", n_models), labels)
  auc_vec       <- setNames(numeric(n_models),        labels)
  best_list     <- setNames(vector("list", n_models), labels)
  legend_labels <- character(n_models)

  for (i in seq_len(n_models)) {
    y_obs  <- models[[i]]$y
    y_pred <- models[[i]]$fitted.values

    roc_df <- roc_for_model(y_obs, y_pred)
    auc    <- trapz_auc(roc_df)
    best   <- roc_df[which.max(roc_df$youden), ]

    roc_list[[i]]  <- roc_df
    auc_vec[i]     <- auc
    best_list[[i]] <- best

    if (ci) {
      cival <- hm_ci(auc, sum(y_obs == 1), sum(y_obs == 0))
      legend_labels[i] <- sprintf(
        "%s: AUC = %.3f (%.3f–%.3f)",
        labels[i], auc, cival["lo"], cival["hi"]
      )
    } else {
      legend_labels[i] <- sprintf("%s: AUC = %.3f", labels[i], auc)
    }
  }

  # ---- Plotting data ---------------------------------------------------------
  plot_df <- do.call(rbind, lapply(seq_len(n_models), function(i) {
    df       <- roc_list[[i]]
    df       <- df[!duplicated(df[, c("fpr", "sensitivity")]), ]
    df$model <- factor(labels[i], levels = labels)
    df
  }))

  best_df <- do.call(rbind, lapply(seq_len(n_models), function(i) {
    b        <- best_list[[i]]
    b$model  <- factor(labels[i], levels = labels)
    b
  }))

  # ---- Subtitle (single model: show core params + AUC) ----------------------
  subtitle <- if (n_models == 1L) {
    coef1    <- models[[1L]]$coefficients
    auc_part <- sub(paste0("^", labels[1L], ": "), "", legend_labels[1L])
    sprintf(
      "E0 = %.2f | Emax = %.2f | EC50 = %.2f\n%s",
      coef1["e0"], coef1["emax"], coef1["ec50"], auc_part
    )
  } else {
    NULL
  }

  # ---- Colour scales ---------------------------------------------------------
  if (is.null(colors)) {
    color_scale <- ggplot2::scale_color_discrete(labels = legend_labels)
    fill_scale  <- ggplot2::scale_fill_discrete(guide = "none")
  } else {
    color_scale <- ggplot2::scale_color_manual(
      values = colors[seq_len(n_models)], labels = legend_labels
    )
    fill_scale  <- ggplot2::scale_fill_manual(
      values = colors[seq_len(n_models)], guide = "none"
    )
  }

  # ---- Build plot ------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = fpr, y = sensitivity, color = model, group = model)
  )

  if (isTRUE(fill)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = sensitivity, fill = model),
      alpha = alpha, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = linewidth) +
    ggplot2::geom_point(
      data = best_df,
      ggplot2::aes(x = fpr, y = sensitivity, color = model),
      size = 3, shape = 16, show.legend = FALSE
    ) +
    color_scale +
    fill_scale +
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
      subtitle = subtitle,
      color    = NULL,
      fill     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      aspect.ratio    = 1,
      legend.position = if (n_models > 1L) "bottom" else "none",
      legend.text     = ggplot2::element_text(size = 11)
    )

  # Best-threshold label: single model only (multi would be too cluttered)
  if (n_models == 1L) {
    p <- p + ggplot2::geom_label(
      data = best_df,
      ggplot2::aes(
        x     = fpr + 0.03,
        y     = sensitivity - 0.05,
        label = sprintf(
          "Thr = %.2f\nSe = %.2f, Sp = %.2f",
          threshold, sensitivity, specificity
        )
      ),
      hjust = 0, size = 3.2, color = "firebrick",
      fill = "white", label.size = 0.3,
      inherit.aes = FALSE
    )
  }

  if (ref_line) {
    p <- p + ggplot2::geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", color = "grey50"
    )
  }

  invisible(structure(
    list(
      plot           = p,
      roc_list       = roc_list,
      roc_df         = roc_list[[1L]],
      auc            = auc_vec,
      best_threshold = best_list
    ),
    class = "roc_logistic_emax"
  ))
}


#' Print a `roc_logistic_emax` object
#'
#' Prints AUC and optimal threshold for each model, then renders the ROC plot.
#'
#' @param x A `roc_logistic_emax` object returned by [roc_logistic_emax()].
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [roc_logistic_emax()]
#'
#' @export
print.roc_logistic_emax <- function(x, ...) {
  if (length(x$auc) == 1L) {
    cat(sprintf("AUC: %.4f\n", x$auc[[1L]]))
    best <- x$best_threshold[[1L]]
    cat(sprintf(
      "Optimal threshold (Youden): %.3f  Se = %.3f  Sp = %.3f\n",
      best$threshold, best$sensitivity, best$specificity
    ))
  } else {
    cat("ROC Curve Comparison\n\n")
    for (nm in names(x$auc)) {
      best <- x$best_threshold[[nm]]
      cat(sprintf(
        "%-30s  AUC = %.4f  |  Youden thr = %.3f  (Se = %.3f, Sp = %.3f)\n",
        nm, x$auc[nm], best$threshold, best$sensitivity, best$specificity
      ))
    }
  }
  print(x$plot)
  invisible(x)
}
