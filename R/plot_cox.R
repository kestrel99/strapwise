#' Forest plot of hazard ratios for a bootstrapped Cox model
#'
#' Plots each coefficient as a point estimate (hazard ratio) with a horizontal
#' bootstrap percentile confidence interval bar, one row per term.
#'
#' @param bc A `boot_cox` object returned by [bootstrap_cox()].
#' @param log_scale Logical; if `TRUE` (default) the x-axis is on the log scale
#'   so that symmetric CI bars reflect log-HR symmetry. Set to `FALSE` for a
#'   linear HR scale.
#' @param point_size Numeric; size of the HR point estimate marker. Default `3`.
#' @param line_size Numeric; line width of the CI bars. Default `0.7`.
#' @param ref_line Logical; draw a vertical reference line at HR = 1. Default
#'   `TRUE`.
#' @param palette Character vector of colours; the first element is used for all
#'   points and CI bars.
#' @param base_theme A complete ggplot2 theme object. Defaults to a
#'   `theme_bw`-based theme.
#'
#' @return A `ggplot` object.
#'
#' @seealso [bootstrap_cox()], [plot_cox_survival()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex + ph.ecog,
#'     data = lung
#'   )
#'   bc <- bootstrap_cox(fit, n_boot = 200, seed = 1)
#'   plot_cox_forest(bc)
#' }
#'
#' @export
plot_cox_forest <- function(
  bc,
  log_scale  = TRUE,
  point_size = 3,
  line_size  = 0.7,
  ref_line   = TRUE,
  palette    = c("#2C7BB6", "#D7191C", "#1A9641", "#F46D43", "#756BB1", "#FDAE61"),
  base_theme = NULL
) {
  if (!inherits(bc, "boot_cox"))
    stop("`bc` must be a boot_cox object returned by bootstrap_cox().")

  alpha  <- 1 - bc$conf_level
  cl_pct <- round(bc$conf_level * 100)

  coef_orig  <- coef(bc$fit)
  ci_lo_log  <- apply(bc$boot_coefs, 2, stats::quantile,
                      probs = alpha / 2,       na.rm = TRUE)
  ci_hi_log  <- apply(bc$boot_coefs, 2, stats::quantile,
                      probs = 1 - alpha / 2,   na.rm = TRUE)

  forest_df <- data.frame(
    term  = factor(names(coef_orig), levels = rev(names(coef_orig))),
    hr    = exp(coef_orig),
    ci_lo = exp(ci_lo_log),
    ci_hi = exp(ci_hi_log),
    stringsAsFactors = FALSE
  )

  use_theme <- if (!is.null(base_theme)) {
    base_theme
  } else {
    ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "none",
        plot.title       = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey45"),
        axis.title       = ggplot2::element_text(size = 10)
      )
  }

  p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = hr, y = term))

  if (ref_line)
    p <- p + ggplot2::geom_vline(
      xintercept = 1, linetype = "dashed", colour = "grey50"
    )

  p <- p +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
      height    = 0.2,
      linewidth = line_size,
      colour    = palette[[1L]]
    ) +
    ggplot2::geom_point(
      size   = point_size,
      shape  = 21,
      fill   = palette[[1L]],
      colour = "grey20"
    ) +
    ggplot2::labs(
      title    = "Hazard Ratio Estimates",
      subtitle = sprintf(
        "%d%% bootstrap percentile CI  (n = %d replicates)",
        cl_pct, bc$n_boot
      ),
      x = "Hazard Ratio",
      y = NULL
    ) +
    use_theme

  if (log_scale)
    p <- p + ggplot2::scale_x_log10()

  p
}


#' Adjusted survival curves for a bootstrapped Cox model
#'
#' Plots model-predicted survival (or event probability, or cumulative hazard)
#' over time for each row of `newdata`, with bootstrap percentile confidence
#' ribbons.
#'
#' Confidence bands are computed by fixing the baseline cumulative hazard at the
#' original-fit estimate and varying only the linear predictor across bootstrap
#' replicates:
#' \deqn{S_b(t \mid x) = \exp\!\bigl(-\hat{H}_0(t)\cdot e^{x^\top\beta_b}\bigr)}
#' where \eqn{\hat{H}_0} is from [survival::basehaz()] and \eqn{\beta_b} is the
#' \eqn{b}-th bootstrap coefficient vector. This correctly reflects uncertainty
#' in the covariate effects; uncertainty in the baseline hazard is not
#' propagated.
#'
#' Stratified Cox models (those using `strata()` in the formula) are not
#' currently supported; an error is raised if strata are detected.
#'
#' @param bc A `boot_cox` object returned by [bootstrap_cox()].
#' @param newdata A data frame where each row defines covariate values for one
#'   predicted curve. Must contain all predictor columns used in the model.
#' @param curve_labels Optional character vector of labels for each row of
#'   `newdata`. Defaults to `"Curve 1"`, `"Curve 2"`, etc.
#' @param times Optional numeric vector of time points at which to evaluate the
#'   curves. `NULL` (default) uses all event times from the original data.
#' @param fun One of `"survival"` (default), `"event"` (1 - S(t)), or
#'   `"cumhaz"` (-log S(t)).
#' @param ribbon_alpha Numeric in `(0, 1)`; CI ribbon transparency. Default
#'   `0.15`.
#' @param line_size Numeric; curve line width. Default `0.9`.
#' @param palette Character vector of colours, one per curve (recycled as
#'   needed).
#' @param base_theme A complete ggplot2 theme object. Defaults to a
#'   `theme_bw`-based theme.
#' @param x_label Character; x-axis label. Default `"Time"`.
#' @param y_label Character; y-axis label. `NULL` (default) is set
#'   automatically from `fun`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [bootstrap_cox()], [plot_cox_forest()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex + ph.ecog,
#'     data = lung
#'   )
#'   bc <- bootstrap_cox(fit, n_boot = 200, seed = 1)
#'
#'   nd <- data.frame(
#'     age    = c(50, 70),
#'     sex    = c(1,  1),
#'     ph.ecog = c(0, 2)
#'   )
#'   plot_cox_survival(bc, newdata = nd,
#'                     curve_labels = c("Age 50, ECOG 0", "Age 70, ECOG 2"))
#' }
#'
#' @export
plot_cox_survival <- function(
  bc,
  newdata,
  curve_labels = NULL,
  times        = NULL,
  fun          = c("survival", "event", "cumhaz"),
  ribbon_alpha = 0.15,
  line_size    = 0.9,
  palette      = c("#2C7BB6", "#D7191C", "#1A9641", "#F46D43", "#756BB1", "#FDAE61"),
  base_theme   = NULL,
  x_label      = "Time",
  y_label      = NULL
) {
  if (!inherits(bc, "boot_cox"))
    stop("`bc` must be a boot_cox object returned by bootstrap_cox().")
  if (!is.data.frame(newdata) || nrow(newdata) < 1L)
    stop("`newdata` must be a non-empty data frame.")

  fun <- match.arg(fun)

  # ---- Baseline cumulative hazard --------------------------------------------
  H0 <- survival::basehaz(bc$fit, centered = FALSE)
  if ("strata" %in% names(H0))
    stop(
      "plot_cox_survival() does not currently support stratified Cox models. ",
      "Use survfit() directly for stratified predictions."
    )

  # Subset to requested times (step-function interpolation)
  if (!is.null(times)) {
    H0_at <- stats::approx(
      x      = H0$time,
      y      = H0$hazard,
      xout   = times,
      method = "constant",  # left-continuous step function
      rule   = 2,           # extrapolate with boundary values
      f      = 0            # left-continuous
    )
    H0 <- data.frame(time = H0_at$x, hazard = H0_at$y)
  }

  t_vals   <- H0$time
  h0_vals  <- H0$hazard
  n_times  <- length(t_vals)
  n_curves <- nrow(newdata)
  n_boot   <- nrow(bc$boot_coefs)
  alpha    <- 1 - bc$conf_level
  cl_pct   <- round(bc$conf_level * 100)

  if (is.null(curve_labels))
    curve_labels <- paste("Curve", seq_len(n_curves))
  if (length(curve_labels) != n_curves)
    stop("`curve_labels` must have one entry per row of `newdata`.")

  # ---- Model matrix for newdata (no intercept, aligned to coef order) -------
  trms <- stats::delete.response(stats::terms(bc$fit))
  mm   <- stats::model.matrix(trms, data = newdata)
  int_col <- which(colnames(mm) == "(Intercept)")
  if (length(int_col) > 0L) mm <- mm[, -int_col, drop = FALSE]
  mm <- mm[, names(coef(bc$fit)), drop = FALSE]

  # ---- Linear predictors (uncentered) ----------------------------------------
  lp_point <- drop(mm %*% coef(bc$fit))                 # n_curves
  lp_boot  <- mm %*% t(bc$boot_coefs)                   # n_curves x n_boot

  # ---- Survival function transformation --------------------------------------
  .transform <- switch(fun,
    survival = function(s) s,
    event    = function(s) 1 - s,
    cumhaz   = function(s) -log(s)
  )

  auto_y_label <- switch(fun,
    survival = "Survival Probability S(t)",
    event    = "Event Probability 1 - S(t)",
    cumhaz   = "Cumulative Hazard H(t)"
  )
  y_label <- if (!is.null(y_label)) y_label else auto_y_label

  # ---- Point estimates and bootstrap CI bands --------------------------------
  # S(t|x) = exp(-H0(t) * exp(lp))
  # outer(exp(lp), h0_vals): n_curves x n_times matrix of H0(t)*exp(lp)
  S_point  <- exp(-outer(exp(lp_point), h0_vals))       # n_curves x n_times

  ci_lo_mat <- matrix(NA_real_, n_curves, n_times)
  ci_hi_mat <- matrix(NA_real_, n_curves, n_times)

  for (i in seq_len(n_curves)) {
    # S_b_i: n_boot x n_times
    S_b_i            <- exp(-outer(exp(lp_boot[i, ]), h0_vals))
    ci_lo_mat[i, ]   <- apply(S_b_i, 2L, stats::quantile,
                               probs = alpha / 2,       na.rm = TRUE)
    ci_hi_mat[i, ]   <- apply(S_b_i, 2L, stats::quantile,
                               probs = 1 - alpha / 2,   na.rm = TRUE)
  }

  # ---- Long-format data frame for ggplot ------------------------------------
  curve_colours <- stats::setNames(rep_len(palette, n_curves), curve_labels)

  plot_df <- do.call(rbind, lapply(seq_len(n_curves), function(i) {
    data.frame(
      time  = t_vals,
      surv  = .transform(S_point[i, ]),
      ci_lo = .transform(ci_lo_mat[i, ]),
      ci_hi = .transform(ci_hi_mat[i, ]),
      curve = curve_labels[[i]],
      stringsAsFactors = FALSE
    )
  }))
  plot_df$curve <- factor(plot_df$curve, levels = curve_labels)

  # For "event" and "cumhaz" the CI bounds are swapped because the
  # transformation is monotone-decreasing (event) or monotone-increasing
  # but the bootstrap quantiles were on the S scale.
  if (fun == "event") {
    lo_tmp        <- plot_df$ci_lo
    plot_df$ci_lo <- plot_df$ci_hi
    plot_df$ci_hi <- lo_tmp
  }

  # ---- Theme -----------------------------------------------------------------
  use_theme <- if (!is.null(base_theme)) {
    base_theme
  } else {
    ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey45"),
        axis.title       = ggplot2::element_text(size = 10),
        legend.title     = ggplot2::element_text(size = 9),
        legend.text      = ggplot2::element_text(size = 8)
      )
  }

  subtitle_txt <- sprintf(
    "%d%% bootstrap CI  (n = %d replicates, fixed baseline hazard)",
    cl_pct, n_boot
  )

  # ---- Plot ------------------------------------------------------------------
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = time, y = surv, colour = curve, fill = curve, group = curve)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      alpha        = ribbon_alpha,
      colour       = NA
    ) +
    ggplot2::geom_step(linewidth = line_size) +
    ggplot2::scale_colour_manual(name = NULL, values = curve_colours) +
    ggplot2::scale_fill_manual(  name = NULL, values = curve_colours) +
    ggplot2::scale_y_continuous(
      name   = y_label,
      limits = if (fun %in% c("survival", "event")) c(0, 1) else NULL,
      breaks = if (fun %in% c("survival", "event")) seq(0, 1, 0.2) else ggplot2::waiver(),
      labels = if (fun %in% c("survival", "event"))
        scales::percent_format(accuracy = 1) else ggplot2::waiver(),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(name = x_label) +
    ggplot2::labs(
      title    = "Adjusted Survival Curves",
      subtitle = subtitle_txt
    ) +
    use_theme
}
