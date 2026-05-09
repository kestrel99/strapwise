#' Response curve plots for a bootstrapped logistic regression
#'
#' Generates one marginal response-curve plot per predictor term in a
#' `boot_logistic` model, using the bootstrap coefficient matrix to construct
#' pointwise confidence intervals. Non-focal predictors are held at their
#' **median** (numeric) or **mode** (factor / character).
#'
#' @details
#' Two plot geometries are used depending on predictor type:
#'
#' - **Continuous predictors**: a probability curve with a shaded bootstrap CI
#'   ribbon. Observed binary outcomes (0 / 1) are overlaid as jittered points
#'   at their true y-value when `raw_data = TRUE`. When `obs_groups` is not
#'   `NULL`, the empirical response rate per quantile stratum is overlaid as
#'   filled squares with Wilson score CI error bars.
#' - **Categorical predictors**: predicted probability as a filled dot with
#'   bootstrap CI error bars per level. Observed 0 / 1 outcomes are overlaid
#'   as jittered points when `raw_data = TRUE`. When `obs_groups` is not
#'   `NULL`, the empirical response rate per category level is also overlaid.
#'
#' **Colour-coding by a grouping variable** (`colour_by`): when a column name
#' is supplied, the raw jitter points are coloured by group level and a legend
#' is added automatically. The observed-proportion squares and their CI bars are
#' **not** affected by `colour_by` by default; set `obs_by_group = TRUE` to
#' also split and colour the proportion overlay by group. The model curve and
#' ribbon are always drawn in a single colour.
#'
#' `colour_by` may reference any column in `bl$data`, including columns that
#' do not appear in the model formula (e.g. a stratification variable).
#'
#' If the formula contains inline transformations such as `log(age)`, plots are
#' produced on the *original* (untransformed) scale; `model.matrix()` applies
#' the transformation internally when computing predictions.
#'
#' @param bl A `boot_logistic` object returned by [bootstrap_logistic()].
#' @param n_grid Integer; number of equally-spaced grid points for continuous
#'   predictor curves. Default `200L`.
#' @param point_size Numeric; size of the predicted-probability marker in
#'   categorical plots. Default `3`.
#' @param line_size Numeric; line width of the probability curve in continuous
#'   plots. Default `0.9`.
#' @param ribbon_alpha Numeric in `(0, 1)`; transparency of the CI ribbon.
#'   Default `0.20`.
#' @param error_width Numeric; cap width of CI error bars for model-predicted
#'   probabilities in categorical plots. Default `0.15`.
#' @param raw_data Logical; overlay observed binary outcomes as jittered points.
#'   Default `TRUE`.
#' @param raw_point_size Numeric; size of raw-data jitter points. Default `1.5`.
#' @param raw_alpha Numeric in `(0, 1)`; transparency of raw-data jitter points.
#'   Default `0.25`.
#' @param jitter_height Numeric; vertical jitter height for raw-data points.
#'   Default `0.02`.
#' @param obs_groups Integer in `[2, 10]` or `NULL`. Number of quantile strata
#'   for the observed-proportion overlay on continuous predictors (2 = halves,
#'   3 = tertiles, 4 = quartiles). Each category level is its own stratum for
#'   factor predictors regardless of this value. `NULL` suppresses the overlay.
#'   Default `4L`.
#' @param obs_colour Colour for the observed-proportion overlay squares and
#'   error bars. Always used when `obs_by_group = FALSE` (the default), and
#'   ignored when `obs_by_group = TRUE`. Default `"#333333"`.
#' @param obs_size Numeric; size of the observed-proportion marker squares.
#'   Default `3`.
#' @param obs_error_width Numeric; cap width of observed-proportion CI error
#'   bars. Default `0.25`.
#' @param colour_by Character string naming a column in `bl$data` by which to
#'   colour-code the raw jitter points. The column need not appear in the model
#'   formula. Coerced to a factor if not already. `NULL` (default) uses a
#'   single colour for all jitter points.
#' @param obs_exclude Character vector of levels of the `colour_by` variable
#'   that should appear in the observed-proportion overlay as a **single**
#'   unquantiled bin (i.e. one square and one set of error bars covering all
#'   observations in that group), rather than being split into `obs_groups`
#'   quantile strata.  This is useful for groups whose predictor values are
#'   constant or near-constant (e.g. a placebo arm with concentration = 0).
#'   Jitter points are never affected.  `obs_exclude` is ignored when
#'   `colour_by = NULL`. Default `NULL` (all groups binned normally).
#' @param obs_by_group Logical; if `TRUE` **and** `colour_by` is set, the
#'   observed-proportion squares and their CI bars are also split and coloured
#'   by group. On continuous plots the per-group squares are dodged horizontally
#'   by `colour_by_dodge` to prevent overlap; on categorical plots
#'   `position_dodge()` is used. Default `FALSE`.
#' @param colour_by_label Character string used as the legend title when
#'   `colour_by` is set. Defaults to the value of `colour_by`.
#' @param group_palette Character vector of colours for the grouping variable
#'   levels, recycled as needed. When `NULL` (default) a built-in qualitative
#'   palette is used.
#' @param colour_by_dodge Positive numeric; horizontal dodge distance for the
#'   per-group observed-proportion markers on **continuous** plots when
#'   `obs_by_group = TRUE`, expressed as a fraction of the predictor's observed
#'   range. Ignored when `obs_by_group = FALSE` or for categorical plots.
#'   Default `0.015`.
#' @param palette Character vector of colours for the model curve and, for
#'   categorical predictors, the predicted-probability fill. Recycled as needed.
#' @param base_theme A complete ggplot2 theme object applied to all plots.
#'   Defaults to a `theme_bw`-based theme.
#' @param free_y Logical; if `FALSE` (default) all plots share a 0-1 y-axis.
#' @param combine Logical; if `TRUE` and **patchwork** is installed, return a
#'   single combined figure. Default `FALSE`.
#'
#' @return A named list of `ggplot` objects (one per predictor), or a
#'   `patchwork` figure when `combine = TRUE`.
#'
#' @seealso [bootstrap_logistic()]
#'
#' @examples
#' df <- data.frame(
#'   y   = rbinom(200, 1, 0.4),
#'   age = rnorm(200, 50, 10),
#'   sex = factor(sample(c("M", "F"), 200, replace = TRUE))
#' )
#' fit <- glm(y ~ age + sex, data = df, family = binomial)
#' bl  <- bootstrap_logistic(fit, n_boot = 200, seed = 1)
#'
#' # Colour jitter points by sex; obs overlay stays single-colour
#' plots <- plot_logistic_curves(bl, colour_by = "sex")
#' plots[["age"]]
#'
#' # Also split the obs proportion overlay by sex
#' plots2 <- plot_logistic_curves(bl, colour_by = "sex", obs_by_group = TRUE)
#' plots2[["age"]]
#'
#' @export
plot_logistic_curves <- function(
  bl,
  n_grid = 200L,
  point_size = 3,
  line_size = 0.9,
  ribbon_alpha = 0.20,
  error_width = 0.15,
  raw_data = TRUE,
  raw_point_size = 1.5,
  raw_alpha = 0.25,
  jitter_height = 0.02,
  obs_groups = 4L,
  obs_colour = "#333333",
  obs_size = 3,
  obs_error_width = 0.25,
  colour_by = NULL,
  obs_exclude = NULL,
  obs_by_group = FALSE,
  colour_by_label = NULL,
  group_palette = NULL,
  colour_by_dodge = 0.015,
  palette = c("#2C7BB6", "#D7191C", "#1A9641", "#F46D43", "#756BB1", "#FDAE61"),
  base_theme = NULL,
  free_y = FALSE,
  combine = FALSE
) {
  # ---- Input validation ------------------------------------------------------
  if (!inherits(bl, "boot_logistic")) {
    stop(
      "`bl` must be a boot_logistic object returned by bootstrap_logistic()."
    )
  }

  if (!is.null(obs_groups)) {
    if (
      !is.numeric(obs_groups) ||
        length(obs_groups) != 1L ||
        obs_groups < 2 ||
        obs_groups > 10 ||
        obs_groups != round(obs_groups)
    ) {
      stop("`obs_groups` must be NULL or a single integer between 2 and 10.")
    }
    obs_groups <- as.integer(obs_groups)
  }

  fit <- bl$fit
  data <- bl$data

  resp_name <- as.character(formula(fit)[[2]])
  pred_terms <- intersect(all.vars(formula(fit)[[3]]), names(data))

  # ---- colour_by setup -------------------------------------------------------
  use_colour_by <- !is.null(colour_by)

  if (use_colour_by) {
    if (!is.character(colour_by) || length(colour_by) != 1L) {
      stop(
        "`colour_by` must be a single character string naming a column in `bl$data`."
      )
    }
    if (!colour_by %in% names(data)) {
      stop(sprintf(
        "`colour_by` column '%s' not found in `bl$data`.",
        colour_by
      ))
    }
    if (colour_by == resp_name) {
      stop("`colour_by` must not be the response variable.")
    }

    grp_vec <- factor(data[[colour_by]])
    grp_lvls <- levels(grp_vec)
    n_grps <- length(grp_lvls)

    default_grp_pal <- c(
      "#E41A1C",
      "#377EB8",
      "#4DAF4A",
      "#984EA3",
      "#FF7F00",
      "#A65628",
      "#F781BF",
      "#999999"
    )
    grp_colours <- stats::setNames(
      if (!is.null(group_palette)) {
        rep_len(group_palette, n_grps)
      } else {
        rep_len(default_grp_pal, n_grps)
      },
      grp_lvls
    )
    legend_title <- if (!is.null(colour_by_label)) {
      colour_by_label
    } else {
      colour_by
    }
    # Centred offsets for horizontal dodge of grouped obs squares
    grp_offsets_raw <- seq_len(n_grps) - (n_grps + 1) / 2
  }

  # Convenience: colour obs overlay by group only when both flags are set
  colour_obs <- use_colour_by && obs_by_group

  # ---- obs_exclude setup -----------------------------------------------------
  # obs_excl_lvls: levels that appear as a single unquantiled bin in the obs
  # overlay.  Jitter is never affected.
  if (!is.null(obs_exclude)) {
    if (!is.character(obs_exclude)) {
      stop("`obs_exclude` must be a character vector of group levels.")
    }
    if (!use_colour_by) {
      warning("`obs_exclude` is ignored when `colour_by = NULL`.")
    }
  }
  obs_excl_lvls <- if (use_colour_by && !is.null(obs_exclude)) {
    obs_exclude
  } else {
    character(0)
  }

  # .obs_split(): compute the obs overlay with per-group binning, where groups
  # in obs_excl_lvls are forced to a single bin (n_groups = 1L triggers the
  # single-stratum fallback in .obs_proportions).  Works for both continuous and
  # categorical predictors; for categoricals the split makes no difference
  # (each level is already its own stratum), but is handled consistently.
  .obs_split <- function(x_sub, y_sub, grp_sub, n_grp, conf_lvl) {
    parts <- list()
    # Non-excluded: normal binning
    keep <- if (length(obs_excl_lvls) > 0L) {
      !as.character(grp_sub) %in% obs_excl_lvls
    } else {
      rep(TRUE, length(x_sub))
    }
    if (any(keep)) {
      parts$normal <- .obs_proportions(
        x_sub[keep],
        y_sub[keep],
        n_groups = n_grp,
        conf_level = conf_lvl,
        group = grp_sub[keep]
      )
    }
    # Excluded: forced single bin per group
    excl <- !keep
    if (any(excl)) {
      parts$single <- .obs_proportions(
        x_sub[excl],
        y_sub[excl],
        n_groups = 1L,
        conf_level = conf_lvl,
        group = grp_sub[excl]
      )
    }
    out <- do.call(rbind, Filter(Negate(is.null), parts))
    if (is.null(out) || nrow(out) == 0L) NULL else out
  }

  cl_pct <- round(bl$conf_level * 100)

  # ---- Default theme ---------------------------------------------------------
  use_theme <- if (!is.null(base_theme)) {
    base_theme
  } else {
    ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 11, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey45"),
        axis.title = ggplot2::element_text(size = 10)
      )
  }

  if (use_colour_by) {
    use_theme <- use_theme +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 9),
        legend.text = ggplot2::element_text(size = 8)
      )
  } else {
    use_theme <- use_theme + ggplot2::theme(legend.position = "none")
  }

  y_limits <- if (free_y) NULL else c(0, 1)
  subtitle_txt <- sprintf(
    "Non-focal predictors at median / mode  |  %d%% bootstrap CI  (n = %d replicates)",
    cl_pct,
    bl$n_boot
  )

  # ---- One plot per predictor term -------------------------------------------
  plots <- vector("list", length(pred_terms))
  names(plots) <- pred_terms

  for (term in pred_terms) {
    col <- data[[term]]
    is_factor <- is.factor(col) || is.character(col)

    pred_df <- .build_pred_frame(
      bl,
      term,
      n_grid = if (is_factor) 1L else n_grid
    )
    pred_res <- .compute_preds(bl, pred_df)

    plot_data <- cbind(pred_df[, term, drop = FALSE], pred_res)
    names(plot_data)[1] <- "x"

    y_label <- paste0("P(", resp_name, " = 1)")

    # =========================================================================
    if (!is_factor) {
      # ---- Continuous predictor -----------------------------------------------

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = prob)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
          fill = palette[1],
          alpha = ribbon_alpha
        ) +
        ggplot2::geom_line(colour = palette[1], linewidth = line_size)

      # -- Raw jitter ----------------------------------------------------------
      if (raw_data) {
        if (use_colour_by) {
          raw_df <- data.frame(
            x = col,
            y = as.numeric(data[[resp_name]]),
            grp = grp_vec,
            stringsAsFactors = FALSE
          )
          raw_df <- raw_df[stats::complete.cases(raw_df), ]
          p <- p +
            ggplot2::geom_jitter(
              data = raw_df,
              ggplot2::aes(x = x, y = y, colour = grp),
              width = 0,
              height = jitter_height,
              size = raw_point_size,
              alpha = raw_alpha,
              inherit.aes = FALSE
            )
        } else {
          raw_df <- data.frame(x = col, y = as.numeric(data[[resp_name]]))
          raw_df <- raw_df[stats::complete.cases(raw_df), ]
          p <- p +
            ggplot2::geom_jitter(
              data = raw_df,
              ggplot2::aes(x = x, y = y),
              width = 0,
              height = jitter_height,
              size = raw_point_size,
              colour = palette[1],
              alpha = raw_alpha,
              inherit.aes = FALSE
            )
        }
      }

      # -- Observed-proportion overlay -----------------------------------------
      if (!is.null(obs_groups)) {
        if (colour_obs) {
          obs_df <- .obs_split(
            col,
            data[[resp_name]],
            grp_vec,
            obs_groups,
            bl$conf_level
          )
          if (!is.null(obs_df)) {
            x_range <- diff(range(col, na.rm = TRUE))
            step <- colour_by_dodge * x_range
            obs_df$grp <- factor(obs_df$grp, levels = grp_lvls)
            obs_df$x <- obs_df$x +
              grp_offsets_raw[as.integer(obs_df$grp)] * step

            p <- p +
              ggplot2::geom_errorbar(
                data = obs_df,
                ggplot2::aes(x = x, ymin = ci_lo, ymax = ci_hi, colour = grp),
                width = obs_error_width,
                linewidth = 0.7,
                inherit.aes = FALSE
              ) +
              ggplot2::geom_point(
                data = obs_df,
                ggplot2::aes(x = x, y = obs, colour = grp, fill = grp),
                size = obs_size,
                shape = 22,
                inherit.aes = FALSE
              )
          }
        } else {
          # Single-colour overlay (default)  -- no group splitting
          obs_df <- .obs_proportions(
            x = col,
            y = data[[resp_name]],
            n_groups = obs_groups,
            conf_level = bl$conf_level
          )
          if (!is.null(obs_df)) {
            p <- p +
              ggplot2::geom_errorbar(
                data = obs_df,
                ggplot2::aes(x = x, ymin = ci_lo, ymax = ci_hi),
                width = obs_error_width,
                linewidth = 0.7,
                colour = obs_colour,
                inherit.aes = FALSE
              ) +
              ggplot2::geom_point(
                data = obs_df,
                ggplot2::aes(x = x, y = obs),
                size = obs_size,
                shape = 22,
                fill = obs_colour,
                colour = obs_colour,
                inherit.aes = FALSE
              )
          }
        }
      }

      # -- Scales and labels ---------------------------------------------------
      p <- p +
        ggplot2::scale_y_continuous(
          name = y_label,
          limits = y_limits,
          breaks = seq(0, 1, 0.2),
          labels = scales::percent_format(accuracy = 1),
          oob = scales::squish
        ) +
        ggplot2::scale_x_continuous(name = term) +
        ggplot2::labs(
          title = paste("Response curve:", term),
          subtitle = subtitle_txt
        ) +
        use_theme

      if (use_colour_by) {
        p <- p +
          ggplot2::scale_colour_manual(
            name = legend_title,
            values = grp_colours
          ) +
          ggplot2::scale_fill_manual(name = legend_title, values = grp_colours)
      }

      # =========================================================================
    } else {
      # ---- Categorical predictor ----------------------------------------------

      lvls <- levels(factor(col))
      fill_colours <- stats::setNames(rep_len(palette, length(lvls)), lvls)
      plot_data$x <- factor(plot_data$x, levels = lvls)

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = prob, fill = x))

      # -- Raw jitter ----------------------------------------------------------
      if (raw_data) {
        if (use_colour_by) {
          raw_df <- data.frame(
            x = factor(col, levels = lvls),
            y = as.numeric(data[[resp_name]]),
            grp = grp_vec,
            stringsAsFactors = FALSE
          )
          raw_df <- raw_df[stats::complete.cases(raw_df), ]
          p <- p +
            ggplot2::geom_jitter(
              data = raw_df,
              ggplot2::aes(x = x, y = y, colour = grp),
              width = 0.15,
              height = jitter_height,
              size = raw_point_size,
              alpha = raw_alpha,
              inherit.aes = FALSE
            )
        } else {
          raw_df <- data.frame(
            x = factor(col, levels = lvls),
            y = as.numeric(data[[resp_name]])
          )
          raw_df <- raw_df[stats::complete.cases(raw_df), ]
          p <- p +
            ggplot2::geom_jitter(
              data = raw_df,
              ggplot2::aes(x = x, y = y, colour = x),
              width = 0.15,
              height = jitter_height,
              size = raw_point_size,
              alpha = raw_alpha,
              inherit.aes = FALSE
            ) +
            ggplot2::scale_colour_manual(values = fill_colours, guide = "none")
        }
      }

      # -- Model predicted probability markers ---------------------------------
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
          width = error_width,
          linewidth = 0.7,
          colour = "grey35"
        ) +
        ggplot2::geom_point(size = point_size, shape = 21, colour = "grey20") +
        ggplot2::scale_fill_manual(values = fill_colours, guide = "none") +
        ggplot2::scale_y_continuous(
          name = y_label,
          limits = y_limits,
          breaks = seq(0, 1, 0.2),
          labels = scales::percent_format(accuracy = 1),
          oob = scales::squish
        ) +
        ggplot2::scale_x_discrete(name = term) +
        ggplot2::labs(
          title = paste("Predicted probability by", term),
          subtitle = subtitle_txt
        ) +
        use_theme

      # -- Observed-proportion overlay -----------------------------------------
      # Squares (shape 22) nudged right of the model circle (shape 21) so the
      # two sources of information are immediately separable at a glance.
      # Error bar cap width is scaled to match the square size: obs_size (mm)
      # divided by a constant that approximates typical mm-per-category spacing,
      # ensuring caps are no wider than the plotting symbol.
      if (!is.null(obs_groups)) {
        cat_eb_width <- obs_size / 30
        nudge <- ggplot2::position_nudge(x = 0.2)
        if (colour_obs) {
          obs_df <- .obs_proportions(
            col,
            data[[resp_name]],
            n_groups = 1L, # ignored for factors
            conf_level = bl$conf_level,
            group = grp_vec
          )
          if (!is.null(obs_df)) {
            obs_df$x <- factor(obs_df$x, levels = lvls)
            obs_df$grp <- factor(obs_df$grp, levels = grp_lvls)
            dodge_nudge <- ggplot2::position_dodge(width = 0.5)
            p <- p +
              ggplot2::geom_errorbar(
                data = obs_df,
                ggplot2::aes(x = x, ymin = ci_lo, ymax = ci_hi, colour = grp),
                width = cat_eb_width,
                linewidth = 0.6,
                position = dodge_nudge,
                inherit.aes = FALSE
              ) +
              ggplot2::geom_point(
                data = obs_df,
                ggplot2::aes(x = x, y = obs, colour = grp, fill = grp),
                size = obs_size,
                shape = 22,
                stroke = 1,
                position = dodge_nudge,
                inherit.aes = FALSE
              ) +
              ggplot2::scale_colour_manual(
                name = legend_title,
                values = grp_colours
              ) +
              ggplot2::scale_fill_manual(
                name = legend_title,
                values = grp_colours
              )
          }
        } else {
          obs_df <- .obs_proportions(
            x = col,
            y = data[[resp_name]],
            n_groups = 1L, # ignored for factors
            conf_level = bl$conf_level
          )
          if (!is.null(obs_df)) {
            obs_df$x <- factor(obs_df$x, levels = lvls)
            p <- p +
              ggplot2::geom_errorbar(
                data = obs_df,
                ggplot2::aes(x = x, ymin = ci_lo, ymax = ci_hi),
                width = cat_eb_width,
                linewidth = 0.6,
                colour = obs_colour,
                position = nudge,
                inherit.aes = FALSE
              ) +
              ggplot2::geom_point(
                data = obs_df,
                ggplot2::aes(x = x, y = obs),
                size = obs_size,
                shape = 22,
                stroke = 1,
                fill = obs_colour,
                colour = "white",
                position = nudge,
                inherit.aes = FALSE
              )
          }
        }
      }

      # When colour_by is active but obs overlay is not grouped, add colour scale
      # (needed for the jitter; fill scale for predicted marker already suppressed)
      if (use_colour_by && !colour_obs) {
        p <- p +
          ggplot2::scale_colour_manual(
            name = legend_title,
            values = grp_colours
          )
      }
    }

    plots[[term]] <- p
  }

  # ---- Optional patchwork combination ----------------------------------------
  if (combine) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(plots))
    }
    message("patchwork is not installed; returning a named list instead.")
  }

  plots
}
