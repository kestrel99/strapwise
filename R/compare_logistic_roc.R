#' Compare logistic regression models with overlaid ROC curves
#'
#' Builds a single ROC plot that overlays the performance of multiple fitted
#' binomial [glm()] models using **ggplot2** and **plotROC**.
#'
#' Models are aligned on the intersection of their fitted observations, using
#' row names from each model frame. This allows comparison of models that were
#' fit on slightly different subsets because of missing predictors, as long as
#' they still share some observations in common and the binary outcome agrees on
#' those shared rows.
#'
#' @param models A list of at least two fitted [glm()] objects with
#'   `family = binomial`.
#' @param model_names Optional character vector of display names for `models`.
#'   Defaults to `names(models)` when present, otherwise `"Model 1"`,
#'   `"Model 2"`, and so on.
#' @param n_cuts Integer; number of cutoff labels drawn along each ROC curve.
#'   Default `0L` for clean overlaid curves.
#' @param labels Logical; whether to draw cutoff labels. Default `FALSE`.
#' @param show_auc Logical; if `TRUE` (default), append AUC values to the legend
#'   labels.
#' @param auc_digits Integer; number of decimal places used when formatting AUC
#'   values in the legend. Default `3L`.
#' @param diagonal Logical; if `TRUE` (default), add a dashed 45-degree
#'   reference line.
#' @param legend_title Character string for the colour legend title. Default
#'   `"Model"`.
#' @param base_theme Optional complete ggplot2 theme object added after
#'   [plotROC::style_roc()]. Default `NULL`.
#'
#' @return A `ggplot` object containing one overlaid ROC curve per model.
#'
#' @seealso [plot_logistic_curves()], [bootstrap_logistic()]
#'
#' @examples
#' df <- data.frame(
#'   y = rbinom(250, 1, 0.5),
#'   age = rnorm(250, 50, 10),
#'   bmi = rnorm(250, 27, 5),
#'   sex = factor(sample(c("F", "M"), 250, replace = TRUE))
#' )
#' df$y <- rbinom(
#'   250, 1,
#'   plogis(-4 + 0.05 * df$age + 0.06 * df$bmi + 0.5 * (df$sex == "M"))
#' )
#'
#' fit_age <- glm(y ~ age, data = df, family = binomial)
#' fit_full <- glm(y ~ age + bmi + sex, data = df, family = binomial)
#'
#' compare_logistic_roc(
#'   list("Age only" = fit_age, "Full model" = fit_full)
#' )
#'
#' @export
compare_logistic_roc <- function(
  models,
  model_names = names(models),
  n_cuts = 0L,
  labels = FALSE,
  show_auc = TRUE,
  auc_digits = 3L,
  diagonal = TRUE,
  legend_title = "Model",
  base_theme = NULL
) {
  if (!is.list(models) || length(models) < 2L) {
    stop("`models` must be a list of at least two fitted binomial glm objects.")
  }

  if (is.null(model_names)) {
    model_names <- paste("Model", seq_along(models))
  }
  if (
    !is.character(model_names) ||
      length(model_names) != length(models) ||
      anyNA(model_names) ||
      any(model_names == "")
  ) {
    stop(
      "`model_names` must be a character vector with one non-empty name per model."
    )
  }
  if (anyDuplicated(model_names)) {
    stop("`model_names` must be unique.")
  }

  if (
    !is.numeric(n_cuts) ||
      length(n_cuts) != 1L ||
      n_cuts < 0 ||
      n_cuts != round(n_cuts)
  ) {
    stop("`n_cuts` must be a single non-negative integer.")
  }
  n_cuts <- as.integer(n_cuts)

  if (!is.logical(labels) || length(labels) != 1L || is.na(labels)) {
    stop("`labels` must be `TRUE` or `FALSE`.")
  }
  if (!is.logical(show_auc) || length(show_auc) != 1L || is.na(show_auc)) {
    stop("`show_auc` must be `TRUE` or `FALSE`.")
  }
  if (
    !is.numeric(auc_digits) ||
      length(auc_digits) != 1L ||
      auc_digits < 0 ||
      auc_digits != round(auc_digits)
  ) {
    stop("`auc_digits` must be a single non-negative integer.")
  }
  auc_digits <- as.integer(auc_digits)

  if (!is.logical(diagonal) || length(diagonal) != 1L || is.na(diagonal)) {
    stop("`diagonal` must be `TRUE` or `FALSE`.")
  }
  if (
    !is.character(legend_title) ||
      length(legend_title) != 1L ||
      is.na(legend_title)
  ) {
    stop("`legend_title` must be a single character string.")
  }

  roc_parts <- lapply(seq_along(models), function(i) {
    .extract_roc_model_data(models[[i]], model_names[[i]])
  })

  common_rows <- Reduce(intersect, lapply(roc_parts, function(x) x$row_id))
  if (length(common_rows) == 0L) {
    stop("`models` do not share any common fitted observations.")
  }

  aligned <- lapply(roc_parts, function(x) {
    x[match(common_rows, x$row_id), , drop = FALSE]
  })

  response_ref <- aligned[[1]]$D
  same_response <- vapply(
    aligned[-1],
    function(x) identical(x$D, response_ref),
    logical(1)
  )
  if (!all(same_response)) {
    stop("All models must have the same binary outcome on shared observations.")
  }

  roc_data <- do.call(rbind, aligned)
  roc_data$model <- factor(roc_data$model, levels = model_names)

  p <- ggplot2::ggplot(
    roc_data,
    ggplot2::aes(d = D, m = M, color = model)
  ) +
    plotROC::geom_roc(n.cuts = n_cuts, labels = labels, show.legend = TRUE) +
    plotROC::style_roc() +
    ggplot2::labs(
      title = "ROC Curve Comparison",
      subtitle = sprintf("Shared observations: %d", length(common_rows)),
      color = legend_title
    )

  if (diagonal) {
    p <- p +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2,
        colour = "grey60"
      )
  }

  if (show_auc) {
    auc_tbl <- plotROC::calc_auc(p)
    auc_col <- intersect(c("model", "name"), names(auc_tbl))[1]
    auc_tbl <- auc_tbl[
      match(model_names, as.character(auc_tbl[[auc_col]])),
      ,
      drop = FALSE
    ]
    auc_labels <- stats::setNames(
      sprintf(
        "%s (AUC = %.*f)",
        model_names,
        auc_digits,
        auc_tbl$AUC
      ),
      model_names
    )
    p <- p +
      ggplot2::scale_colour_discrete(
        labels = function(x) unname(auc_labels[x])
      )
  }

  if (!is.null(base_theme)) {
    p <- p + base_theme
  }

  p
}


.extract_roc_model_data <- function(fit, model_name) {
  if (!inherits(fit, "glm")) {
    stop("Each element of `models` must be a fitted glm object.")
  }
  if (fit$family$family != "binomial") {
    stop("Each element of `models` must use family = binomial.")
  }

  mf <- stats::model.frame(fit)
  y <- stats::model.response(mf)
  if (is.matrix(y)) {
    stop(
      "ROC comparison requires a binary 0/1 response, not an aggregated binomial response."
    )
  }

  row_id <- rownames(mf)
  if (is.null(row_id)) {
    row_id <- names(stats::fitted(fit))
  }
  if (is.null(row_id)) {
    row_id <- as.character(seq_len(nrow(mf)))
  }

  data.frame(
    row_id = row_id,
    D = as.numeric(plotROC::verify_d(y)),
    M = as.numeric(stats::fitted(fit)),
    model = model_name,
    stringsAsFactors = FALSE
  )
}
