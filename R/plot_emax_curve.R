#' Response curve plot for a bootstrapped logistic Emax model
#'
#' Generates a response probability curve with a bootstrap percentile CI ribbon
#' and optional overlays of observed binary outcomes and binned observed
#' proportions from a `boot_emax` object.
#'
#' Other covariates (linear, Emax modifier, and EC50 modifier) are fixed at
#' their training-data column means for the fitted curve unless custom matrices
#' are supplied via `newlinear`, `newemax`, or `newec50`. When such covariates
#' are present the subtitle reports what values were used.
#'
#' Bootstrap CIs are percentile intervals: each replicate in `be$boot_theta`
#' is used to compute a response curve, and the `alpha/2` and
#' `1 - alpha/2` quantiles across replicates form the ribbon bounds.
#'
#' @param be A `boot_emax` object returned by [bootstrap_emax()].
#' @param n_grid Integer; number of equally-spaced grid points for the fitted
#'   curve. Default `200L`.
#' @param x_range Numeric vector `c(lo, hi)` defining the x-axis range.
#'   Defaults to the observed range in the training data.
#' @param newlinear Optional numeric matrix with `n_grid` rows supplying linear
#'   covariate values for the prediction grid. Must have the same number of
#'   columns as the training `linear_covs`. `NULL` (default) uses training
#'   column means.
#' @param newemax Optional numeric matrix of Emax covariate values for the
#'   prediction grid. Same convention as `newlinear`.
#' @param newec50 Optional numeric matrix of EC50 covariate values for the
#'   prediction grid. Same convention as `newlinear`.
#' @param line_size Numeric; probability curve line width. Default `0.9`.
#' @param ribbon_alpha Numeric in `(0, 1)`; CI ribbon transparency. Default
#'   `0.20`.
#' @param raw_data Logical; overlay observed 0/1 outcomes as jittered points.
#'   Default `TRUE`.
#' @param raw_point_size Numeric; jitter point size. Default `1.5`.
#' @param raw_alpha Numeric in `(0, 1)`; jitter point transparency. Default
#'   `0.25`.
#' @param jitter_height Numeric; vertical jitter height. Default `0.02`.
#' @param obs_groups Integer in `[2, 10]` or `NULL`; number of quantile bins
#'   for the observed-proportion overlay (squares with Wilson score CI bars).
#'   `NULL` suppresses the overlay. Default `4L`.
#' @param obs_colour Colour for the observed-proportion markers and error bars.
#'   Default `"#333333"`.
#' @param obs_size Numeric; size of the observed-proportion squares. Default
#'   `3`.
#' @param obs_error_width Numeric; cap width of observed-proportion CI error
#'   bars. Default `0.25`.
#' @param palette Character vector of colours; the first element is used for
#'   the curve, ribbon, and raw jitter points.
#' @param base_theme A complete ggplot2 theme object applied to the plot.
#'   Defaults to a `theme_bw`-based theme.
#' @param free_y Logical; if `FALSE` (default) the y-axis spans 0-1.
#' @param x_label Character; x-axis label. Default `"x"`.
#' @param y_label Character; y-axis label. Default `"P(Response = 1)"`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [bootstrap_emax()], [fit_logistic_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#' fit  <- fit_logistic_emax(y, dose)
#' be   <- bootstrap_emax(fit, n_boot = 200, seed = 1)
#'
#' plot_emax_curve(be)
#'
#' @export
plot_emax_curve <- function(
  be,
  n_grid          = 200L,
  x_range         = NULL,
  newlinear       = NULL,
  newemax         = NULL,
  newec50         = NULL,
  line_size       = 0.9,
  ribbon_alpha    = 0.20,
  raw_data        = TRUE,
  raw_point_size  = 1.5,
  raw_alpha       = 0.25,
  jitter_height   = 0.02,
  obs_groups      = 4L,
  obs_colour      = "#333333",
  obs_size        = 3,
  obs_error_width = 0.25,
  palette         = c("#2C7BB6", "#D7191C", "#1A9641", "#F46D43", "#756BB1", "#FDAE61"),
  base_theme      = NULL,
  free_y          = FALSE,
  x_label         = "x",
  y_label         = "P(Response = 1)"
) {
  if (!inherits(be, "boot_emax"))
    stop("`be` must be a boot_emax object returned by bootstrap_emax().")

  if (!is.null(obs_groups)) {
    if (
      !is.numeric(obs_groups) || length(obs_groups) != 1L ||
      obs_groups < 2 || obs_groups > 10 || obs_groups != round(obs_groups)
    )
      stop("`obs_groups` must be NULL or a single integer between 2 and 10.")
    obs_groups <- as.integer(obs_groups)
  }

  fit    <- be$fit
  idx    <- fit$idx
  alpha  <- 1 - be$conf_level
  cl_pct <- round(be$conf_level * 100)

  # ---- x grid ----------------------------------------------------------------
  lo     <- if (!is.null(x_range)) x_range[[1L]] else min(fit$x, na.rm = TRUE)
  hi     <- if (!is.null(x_range)) x_range[[2L]] else max(fit$x, na.rm = TRUE)
  x_grid <- seq(lo, hi, length.out = n_grid)

  # ---- Covariate matrices at prediction points ------------------------------
  make_grid_cov <- function(train_z, user_z, label) {
    p <- ncol(train_z)
    if (p == 0L) return(matrix(0, nrow = n_grid, ncol = 0L))
    if (!is.null(user_z)) {
      m <- as.matrix(user_z)
      if (nrow(m) != n_grid) stop(label, " must have ", n_grid, " rows.")
      if (ncol(m) != p)      stop(label, " must have ", p, " column(s).")
      return(m)
    }
    matrix(colMeans(train_z), nrow = n_grid, ncol = p, byrow = TRUE)
  }

  Z_lin_p  <- make_grid_cov(fit$Z_lin,  newlinear, "newlinear")
  Z_emax_p <- make_grid_cov(fit$Z_emax, newemax,   "newemax")
  Z_ec50_p <- make_grid_cov(fit$Z_ec50, newec50,   "newec50")

  p_lin  <- ncol(Z_lin_p)
  p_emax <- ncol(Z_emax_p)
  p_ec50 <- ncol(Z_ec50_p)

  # ---- Prediction closure: one theta vector -> probability vector over x_grid
  .pred <- function(th) {
    emax_i <- th[idx$emax0] +
      if (p_emax > 0L) drop(Z_emax_p %*% th[idx$emax_cov]) else 0
    ec50_i <- exp(th[idx$log_ec50] +
      if (p_ec50 > 0L) drop(Z_ec50_p %*% th[idx$ec50_cov]) else 0)
    eta <- th[idx$e0] + emax_i * x_grid / (ec50_i + x_grid) +
      if (p_lin > 0L) drop(Z_lin_p %*% th[idx$lin]) else 0
    as.vector(stats::plogis(eta))
  }

  # ---- Point estimate and bootstrap CI band ----------------------------------
  prob     <- .pred(fit$theta_hat)
  # apply over rows of boot_theta (each row is one theta vector)
  # result is n_grid x n_boot
  boot_mat <- apply(be$boot_theta, 1L, .pred)
  ci_lo    <- apply(boot_mat, 1L, stats::quantile, probs = alpha / 2,       na.rm = TRUE)
  ci_hi    <- apply(boot_mat, 1L, stats::quantile, probs = 1 - alpha / 2,   na.rm = TRUE)

  curve_df <- data.frame(x = x_grid, prob = prob, ci_lo = ci_lo, ci_hi = ci_hi)

  # ---- Theme -----------------------------------------------------------------
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

  y_limits <- if (free_y) NULL else c(0, 1)

  has_covs   <- p_lin > 0L || p_emax > 0L || p_ec50 > 0L
  covs_fixed <- any(!is.null(newlinear), !is.null(newemax), !is.null(newec50))
  subtitle_txt <- sprintf(
    "%s%d%% bootstrap CI  (n = %d replicates)",
    if (has_covs)
      paste0(
        "Other covariates at ",
        if (covs_fixed) "supplied values" else "training column means",
        "  |  "
      )
    else "",
    cl_pct, be$n_boot
  )

  # ---- Base plot -------------------------------------------------------------
  p <- ggplot2::ggplot(curve_df, ggplot2::aes(x = x, y = prob)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      fill  = palette[[1L]],
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(colour = palette[[1L]], linewidth = line_size)

  # ---- Raw jitter points -----------------------------------------------------
  if (raw_data) {
    raw_df <- data.frame(x = fit$x, y = as.numeric(fit$y))
    raw_df <- raw_df[stats::complete.cases(raw_df), ]
    p <- p +
      ggplot2::geom_jitter(
        data        = raw_df,
        ggplot2::aes(x = x, y = y),
        width       = 0,
        height      = jitter_height,
        size        = raw_point_size,
        colour      = palette[[1L]],
        alpha       = raw_alpha,
        inherit.aes = FALSE
      )
  }

  # ---- Observed-proportion overlay -------------------------------------------
  if (!is.null(obs_groups)) {
    obs_df <- .obs_proportions(
      x          = fit$x,
      y          = fit$y,
      n_groups   = obs_groups,
      conf_level = be$conf_level
    )
    if (!is.null(obs_df)) {
      p <- p +
        ggplot2::geom_errorbar(
          data        = obs_df,
          ggplot2::aes(x = x, ymin = ci_lo, ymax = ci_hi),
          width       = obs_error_width,
          linewidth   = 0.7,
          colour      = obs_colour,
          inherit.aes = FALSE
        ) +
        ggplot2::geom_point(
          data        = obs_df,
          ggplot2::aes(x = x, y = obs),
          size        = obs_size,
          shape       = 22,
          fill        = obs_colour,
          colour      = obs_colour,
          inherit.aes = FALSE
        )
    }
  }

  # ---- Scales and labels -----------------------------------------------------
  p +
    ggplot2::scale_y_continuous(
      name   = y_label,
      limits = y_limits,
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1),
      oob    = scales::squish
    ) +
    ggplot2::scale_x_continuous(name = x_label) +
    ggplot2::labs(
      title    = "Logistic Emax Response Curve",
      subtitle = subtitle_txt
    ) +
    use_theme
}
