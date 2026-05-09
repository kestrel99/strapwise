#' Overlaid ROC curves for one or more fitted models
#'
#' Computes ROC curves for one or more fitted models using the same trapezoidal
#' AUC and Youden-threshold logic as [roc_logistic_emax()], overlays them on a
#' single panel, and returns ROC data, AUC values, and the ggplot object.
#'
#' Any model class is supported provided `predict(object, type = "response")`
#' returns a numeric vector of fitted probabilities the same length as
#' `outcome`.
#'
#' @param models A named list of one or more fitted model objects.
#' @param outcome Numeric binary (0/1) vector of observed outcomes, the same
#'   length as the fitted values returned by each model.
#' @param labels Optional character vector of display labels, one per model.
#'   Defaults to `names(models)`; unnamed models are labelled `"Model 1"`,
#'   `"Model 2"`, etc.
#' @param thresholds Numeric vector of classification thresholds. When `NULL`
#'   (default), 501 evenly spaced values from 1 to 0 are used.
#' @param ci Logical; if `TRUE` (default) a 95% Hanley-McNeil confidence
#'   interval is appended to each legend entry.
#' @param ref_line Logical; if `TRUE` (default) a dashed diagonal chance line
#'   is drawn.
#' @param fill Logical; whether to shade the area under the ROC curve. When
#'   `NULL` (default), shading is applied for a single model and suppressed for
#'   multiple models.
#' @param alpha Opacity of the AUC shading (single-model case). Default `0.2`.
#' @param colors Optional character vector of line colours, one per model.
#'   Defaults to the ggplot2 discrete colour scale.
#' @param linewidth Line width. Default `1.1`.
#' @param title Plot title. Default `"ROC Curves"`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return An S3 object of class `"roc_logistic"`, a named list with:
#'   \describe{
#'     \item{`plot`}{A `ggplot` object.}
#'     \item{`roc_list`}{Named list of ROC data frames (one per model), each
#'       with columns `threshold`, `sensitivity`, `specificity`, `fpr`, and
#'       `youden`.}
#'     \item{`auc`}{Named numeric vector of AUC values.}
#'     \item{`best_threshold`}{Named list of single-row data frames for the
#'       Youden-optimal threshold of each model.}
#'   }
#'
#' @seealso [roc_logistic_emax()], [compare_logistic_roc()]
#'
#' @examples
#' set.seed(1)
#' n   <- 300
#' age <- rnorm(n, 55, 10)
#' bmi <- rnorm(n, 27, 5)
#' y   <- rbinom(n, 1, plogis(-4 + 0.05 * age + 0.06 * bmi))
#' df  <- data.frame(y = y, age = age, bmi = bmi)
#'
#' m1 <- glm(y ~ age,       data = df, family = binomial)
#' m2 <- glm(y ~ age + bmi, data = df, family = binomial)
#'
#' roc <- roc_logistic(list(Age = m1, Full = m2), outcome = y)
#' print(roc)
#'
#' @export
roc_logistic <- function(
  models,
  outcome,
  labels    = NULL,
  thresholds = NULL,
  ci        = TRUE,
  ref_line  = TRUE,
  fill      = NULL,
  alpha     = 0.2,
  colors    = NULL,
  linewidth = 1.1,
  title     = "ROC Curves",
  ...
) {
  if (!is.list(models)) models <- list(models)
  n_models <- length(models)
  if (n_models == 0L) stop("`models` must contain at least one model.")

  # Labels
  if (is.null(labels)) labels <- names(models)
  if (is.null(labels)) labels <- character(n_models)
  unnamed <- nchar(labels) == 0 | is.na(labels)
  labels[unnamed] <- paste0("Model ", which(unnamed))
  if (length(labels) != n_models)
    stop("`labels` must have the same length as `models`.")
  if (anyDuplicated(labels))
    stop("`labels` must be unique.")

  outcome <- as.numeric(outcome)
  if (!all(outcome %in% c(0, 1, NA)))
    stop("`outcome` must be binary (0/1).")

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
    q1 <- auc / (2 - auc)
    q2 <- 2 * auc^2 / (1 + auc)
    se <- sqrt(
      (auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) /
      (n1 * n0)
    )
    z  <- stats::qnorm(0.975)
    c(lo = pmax(0, auc - z * se), hi = pmin(1, auc + z * se))
  }

  # ---- Compute ROC per model -------------------------------------------------
  roc_list      <- setNames(vector("list", n_models), labels)
  auc_vec       <- setNames(numeric(n_models), labels)
  best_list     <- setNames(vector("list", n_models), labels)
  legend_labels <- character(n_models)

  for (i in seq_len(n_models)) {
    y_pred <- tryCatch(
      predict(models[[i]], type = "response"),
      error = function(e) stop(
        "Could not get fitted probabilities from model '", labels[i],
        "': ", conditionMessage(e)
      )
    )
    y_obs <- outcome[seq_along(y_pred)]
    ok    <- !is.na(y_obs) & !is.na(y_pred)
    if (sum(ok) < 2L)
      stop("Model '", labels[i], "' has fewer than 2 complete observations.")
    y_obs  <- y_obs[ok]
    y_pred <- y_pred[ok]
    if (length(y_pred) != sum(ok))
      stop("Length of `outcome` (", length(outcome), ") does not match the ",
           "number of fitted values for model '", labels[i], "' (",
           length(y_pred), ").")

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

  # ---- Build long-format plotting data ---------------------------------------
  plot_df <- do.call(rbind, lapply(seq_len(n_models), function(i) {
    df <- roc_list[[i]]
    df <- df[!duplicated(df[, c("fpr", "sensitivity")]), ]
    df$model <- factor(labels[i], levels = labels)
    df
  }))

  best_df <- do.call(rbind, lapply(seq_len(n_models), function(i) {
    b <- best_list[[i]]
    b$model <- factor(labels[i], levels = labels)
    b
  }))

  # ---- Colour scales ---------------------------------------------------------
  if (is.null(colors)) {
    color_scale <- ggplot2::scale_color_discrete(labels = legend_labels)
    fill_scale  <- ggplot2::scale_fill_discrete(guide = "none")
  } else {
    color_scale <- ggplot2::scale_color_manual(
      values = colors[seq_len(n_models)], labels = legend_labels
    )
    fill_scale <- ggplot2::scale_fill_manual(
      values = colors[seq_len(n_models)], guide = "none"
    )
  }

  # ---- Plot ------------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = fpr, y = sensitivity, color = model, group = model)
  )

  if (fill) {
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
    { if (ref_line)
        ggplot2::geom_abline(
          slope = 1, intercept = 0, linetype = "dashed", color = "grey50"
        )
    } +
    color_scale +
    fill_scale +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(
      x     = "1 - Specificity (FPR)",
      y     = "Sensitivity (TPR)",
      title = title,
      color = NULL,
      fill  = NULL
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      aspect.ratio    = 1,
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 11)
    )

  invisible(structure(
    list(plot = p, roc_list = roc_list, auc = auc_vec, best_threshold = best_list),
    class = "roc_logistic"
  ))
}


#' Print a `roc_logistic` object
#'
#' Prints AUC and optimal threshold for each model, then renders the ROC plot.
#'
#' @param x A `roc_logistic` object returned by [roc_logistic()].
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [roc_logistic()]
#'
#' @export
print.roc_logistic <- function(x, ...) {
  cat("ROC Curve Comparison\n\n")
  for (nm in names(x$auc)) {
    best <- x$best_threshold[[nm]]
    cat(sprintf(
      "%-30s  AUC = %.4f  |  Youden thr = %.3f  (Se = %.3f, Sp = %.3f)\n",
      nm, x$auc[nm],
      best$threshold, best$sensitivity, best$specificity
    ))
  }
  print(x$plot)
  invisible(x)
}
