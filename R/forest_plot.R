# Internal formatting helpers --------------------------------------------------

# Format to `digits` significant figures; trimws() removes the leading space
# that formatC inserts when right-aligning numeric output.
.fp_fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "-",
         trimws(formatC(signif(x, digits), format = "g", digits = digits)))
}

# Format p-values: "" for NA, "<0.001" for very small values.
.fp_fmt_p <- function(x) {
  ifelse(is.na(x), "",
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

# For an interaction term, find all model coefficients that belong to it.
#
# Matches by splitting both the term label and each coefficient name on ":"
# and checking that every component of the term label is a prefix of some
# component of the coefficient name.  This handles both "age:sexMale"
# (numeric first) and "sexMale:age" (factor first) layouts.
.fp_int_coefs <- function(term, coef_names) {
  vars <- strsplit(term, ":", fixed = TRUE)[[1L]]
  coef_names[vapply(coef_names, function(cn) {
    parts <- strsplit(cn, ":", fixed = TRUE)[[1L]]
    length(parts) == length(vars) &&
      all(vapply(vars, function(v) any(startsWith(parts, v)), logical(1L)))
  }, logical(1L))]
}

# From an interaction coefficient name (e.g. "age:sexMale"), derive a short
# display label by stripping the matched variable name from each component,
# returning only the non-trivial suffixes (e.g. "Male").
.fp_int_label <- function(cn, vars) {
  parts <- strsplit(cn, ":", fixed = TRUE)[[1L]]
  stripped <- vapply(parts, function(p) {
    match_var <- vars[startsWith(p, vars)]
    if (length(match_var)) sub(paste0("^", match_var[[1L]]), "", p) else p
  }, character(1L))
  paste(stripped[nzchar(stripped)], collapse = ":")
}


# fp_data ----------------------------------------------------------------------

#' Build a tidy data frame for a forest plot
#'
#' Prepares a structured data frame from a fitted binomial `glm` or `coxph`
#' object suitable for passing to [fp_plot()].  The model type is
#' auto-detected from the object class.
#'
#' For categorical predictors one row is produced per factor level, with the
#' first level as the reference category.  The p-value on the header row is
#' derived from the Wald statistic of the first non-reference contrast
#' (appropriate for binary predictors; compute a likelihood-ratio p-value
#' separately for predictors with more than two levels).
#'
#' For continuous predictors each row represents an estimated OR/HR relative
#' to the predictor's median in `data`.
#'
#' Interaction terms (those containing `":"`) are detected automatically.  A
#' continuous x continuous interaction is shown as a single row.  An
#' interaction involving a factor is shown as a header row with one child row
#' per factor level present in the model.
#'
#' @param model   A fitted `glm` (binomial family) or `coxph` object.
#' @param data    Data frame used to fit `model`.
#' @param labels  Named character vector mapping raw variable names to
#'   publication-ready display labels.  Unmatched names fall back to the
#'   raw variable name.  Interaction term names (e.g. `"age:bmi"`) can also
#'   be supplied as keys.
#' @param log_vars Character vector of variable names stored on the log scale.
#'   Their percentile display values are back-transformed with `exp()` so
#'   readers see the original units.
#' @param percs   Numeric vector of percentiles shown for continuous
#'   predictors (default `c(0.05, 0.25, 0.75, 0.95)`).
#' @param ci_level Confidence level for the CI text column (default `0.95`).
#' @param digits  Significant figures used when formatting percentile display
#'   values and the estimate/CI text column (default `3`).
#' @param ci_width Integer; number of spaces in the CI-panel spacer column.
#'   Increase to widen the CI drawing area in [fp_plot()] (default `34`).
#'
#' @return A data frame with one row per displayed entry containing:
#'   `Predictor`, `est`, `lci`, `uci`, `pval`, `is_header`, `is_reference`,
#'   `group`, the model-specific estimate column (`OR (95% CI)` or
#'   `HR (95% CI)`), `P-value`, and `" "` (spacer column).
#'
#' @seealso [fp_plot()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   lung$sex <- factor(lung$sex, levels = 1:2, labels = c("Male", "Female"))
#'   lung_cc <- lung[
#'     complete.cases(lung[, c("time", "status", "age", "sex")]), ]
#'
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex,
#'     data = lung_cc
#'   )
#'   pd <- fp_data(fit, lung_cc,
#'                 labels = c(age = "Age (years)", sex = "Sex"))
#'   head(pd[, c("Predictor", "HR (95% CI)", "P-value")])
#' }
#'
#' @export
fp_data <- function(
  model,
  data,
  labels   = NULL,
  log_vars = NULL,
  percs    = c(0.05, 0.25, 0.75, 0.95),
  ci_level = 0.95,
  digits   = 3L,
  ci_width = 34L
) {
  is_cox    <- inherits(model, "coxph")
  est_label <- if (is_cox) "HR (95% CI)" else "OR (95% CI)"

  coefs      <- stats::coef(model)
  vcov_mat   <- stats::vcov(model)
  coef_names <- names(coefs)
  z_crit     <- stats::qnorm(1 - (1 - ci_level) / 2)
  terms_     <- attr(stats::terms(model), "term.labels")

  # -- Input validation ------------------------------------------------------
  if (!is.null(labels)) {
    bad_lbl <- setdiff(names(labels), terms_)
    if (length(bad_lbl))
      warning("fp_data(): label key(s) not found in model terms: ",
              paste(bad_lbl, collapse = ", "),
              call. = FALSE)
  }

  if (!is.null(log_vars)) {
    bad_lv <- setdiff(log_vars, terms_)
    if (length(bad_lv))
      warning("fp_data(): log_vars entry/entries not found in model terms: ",
              paste(bad_lv, collapse = ", "),
              call. = FALSE)
    cat_lv <- intersect(log_vars, terms_)
    cat_lv <- cat_lv[vapply(cat_lv, function(v) {
      col <- data[[v]]
      !is.null(col) && (is.factor(col) || is.character(col))
    }, logical(1L))]
    if (length(cat_lv))
      warning("fp_data(): log_vars contains categorical predictor(s) ",
              "(back-transform has no effect): ",
              paste(cat_lv, collapse = ", "),
              call. = FALSE)
  }

  .lbl  <- function(v) {
    if (!is.null(labels) && v %in% names(labels)) labels[[v]] else v
  }
  .ilog <- function(v) !is.null(log_vars) && v %in% log_vars

  rows <- list()
  grp  <- 0L

  for (term in terms_) {
    grp <- grp + 1L
    is_interaction <- grepl(":", term, fixed = TRUE)
    col    <- if (!is_interaction) data[[term]] else NULL
    is_cat <- !is_interaction && (is.factor(col) || is.character(col))

    if (is_interaction) {
      # -- Interaction --------------------------------------------------------
      int_coefs <- .fp_int_coefs(term, coef_names)
      if (length(int_coefs) == 0L) next
      vars_i <- strsplit(term, ":", fixed = TRUE)[[1L]]

      # A factor interaction has a level suffix after stripping variable names
      # (e.g. "age:sexM" -> suffix "M").  Use this to distinguish from a pure
      # numeric interaction (e.g. "age:bmi" -> suffix ""), which has no levels.
      has_factor_lvl <- any(vapply(int_coefs, function(cn) {
        nzchar(.fp_int_label(cn, vars_i))
      }, logical(1L)))

      if (!has_factor_lvl) {
        # Single row, no header (continuous x continuous)
        b <- coefs[[int_coefs]]
        s <- sqrt(vcov_mat[[int_coefs, int_coefs]])
        rows <- c(rows, list(data.frame(
          Predictor    = .lbl(term),
          est          = exp(b),
          lci          = exp(b - z_crit * s),
          uci          = exp(b + z_crit * s),
          pval         = 2 * stats::pnorm(-abs(b / s)),
          is_header    = FALSE,
          is_reference = FALSE,
          group        = grp,
          stringsAsFactors = FALSE
        )))
      } else {
        # Header + per-level rows (factor involved, one or more non-ref levels)
        b1 <- coefs[[int_coefs[[1L]]]]
        s1 <- sqrt(vcov_mat[[int_coefs[[1L]], int_coefs[[1L]]]])
        rows <- c(rows, list(data.frame(
          Predictor    = .lbl(term),
          est          = NA_real_,
          lci          = NA_real_,
          uci          = NA_real_,
          pval         = 2 * stats::pnorm(-abs(b1 / s1)),
          is_header    = TRUE,
          is_reference = FALSE,
          group        = grp,
          stringsAsFactors = FALSE
        )))
        for (cn in int_coefs) {
          b   <- coefs[[cn]]
          s   <- sqrt(vcov_mat[[cn, cn]])
          lbl <- .fp_int_label(cn, vars_i)
          rows <- c(rows, list(data.frame(
            Predictor    = paste0("  ", if (nzchar(lbl)) lbl else cn),
            est          = exp(b),
            lci          = exp(b - z_crit * s),
            uci          = exp(b + z_crit * s),
            pval         = NA_real_,
            is_header    = FALSE,
            is_reference = FALSE,
            group        = grp,
            stringsAsFactors = FALSE
          )))
        }
      }

    } else if (is_cat) {
      # -- Categorical ------------------------------------------------------
      levs   <- levels(factor(col))
      cn_hdr <- paste0(term, levs[[2L]])
      b_hdr  <- coefs[[cn_hdr]]
      se_hdr <- sqrt(vcov_mat[[cn_hdr, cn_hdr]])
      p_hdr  <- 2 * stats::pnorm(-abs(b_hdr / se_hdr))

      rows <- c(rows, list(data.frame(
        Predictor    = .lbl(term),
        est          = NA_real_,
        lci          = NA_real_,
        uci          = NA_real_,
        pval         = p_hdr,
        is_header    = TRUE,
        is_reference = FALSE,
        group        = grp,
        stringsAsFactors = FALSE
      )))
      rows <- c(rows, list(data.frame(
        Predictor    = paste0("  ", levs[[1L]], " (reference)"),
        est          = NA_real_,
        lci          = NA_real_,
        uci          = NA_real_,
        pval         = NA_real_,
        is_header    = FALSE,
        is_reference = TRUE,
        group        = grp,
        stringsAsFactors = FALSE
      )))
      for (i in seq(2L, length(levs))) {
        cn <- paste0(term, levs[[i]])
        if (!cn %in% coef_names) next
        b  <- coefs[[cn]]
        s  <- sqrt(vcov_mat[[cn, cn]])
        rows <- c(rows, list(data.frame(
          Predictor    = paste0("  ", levs[[i]]),
          est          = exp(b),
          lci          = exp(b - z_crit * s),
          uci          = exp(b + z_crit * s),
          pval         = NA_real_,
          is_header    = FALSE,
          is_reference = FALSE,
          group        = grp,
          stringsAsFactors = FALSE
        )))
      }

    } else {
      # -- Continuous -------------------------------------------------------
      if (!term %in% coef_names) next
      b <- coefs[[term]]
      s <- sqrt(vcov_mat[[term, term]])
      p <- 2 * stats::pnorm(-abs(b / s))

      rows <- c(rows, list(data.frame(
        Predictor    = .lbl(term),
        est          = NA_real_,
        lci          = NA_real_,
        uci          = NA_real_,
        pval         = p,
        is_header    = TRUE,
        is_reference = FALSE,
        group        = grp,
        stringsAsFactors = FALSE
      )))

      med <- stats::median(col, na.rm = TRUE)
      for (perc in percs) {
        pv      <- stats::quantile(col, probs = perc, na.rm = TRUE)
        d       <- pv - med
        bch     <- b * d
        sch     <- s * abs(d)
        dv      <- if (.ilog(term)) exp(pv) else pv
        pct_lbl <- sprintf("%dth percentile", as.integer(perc * 100))
        rows <- c(rows, list(data.frame(
          Predictor    = paste0("  ", pct_lbl, " (", .fp_fmt(dv, digits), ")"),
          est          = exp(bch),
          lci          = exp(bch - z_crit * sch),
          uci          = exp(bch + z_crit * sch),
          pval         = NA_real_,
          is_header    = FALSE,
          is_reference = abs(d) < 1e-10,
          group        = grp,
          stringsAsFactors = FALSE
        )))
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  ci_str <- ifelse(
    is.na(out$est),
    "",
    paste0(
      .fp_fmt(out$est, digits), " (",
      .fp_fmt(out$lci, digits), ", ",
      .fp_fmt(out$uci, digits), ")"
    )
  )

  out[[est_label]] <- ci_str
  out[["P-value"]] <- .fp_fmt_p(out$pval)
  out[[" "]]       <- strrep(" ", ci_width)

  out
}


# fp_plot ----------------------------------------------------------------------

#' Build and render a forest plot from a fitted model
#'
#' Calls [fp_data()] to prepare the plot data frame and then renders a
#' `forestploter` forest plot in one step.  All [fp_data()] arguments are
#' forwarded directly; call [fp_data()] separately first only if you need to
#' inspect or modify the intermediate data frame before plotting.
#'
#' The model type (logistic or Cox) is auto-detected from the object class.
#'
#' The four display columns are:
#' 1. `Predictor`  -- row labels
#' 2. `OR (95% CI)` or `HR (95% CI)`  -- formatted estimate string
#' 3. `" "`  -- spacer column where the CI glyph is drawn
#' 4. `P-value`
#'
#' Header rows (one per predictor) receive a shaded background and bold label.
#' P-values below 0.05 are coloured red.
#'
#' @param model   A fitted `glm` (binomial family) or `coxph` object.
#' @param data    Data frame used to fit `model`.
#' @param labels  Named character vector; display labels for predictor names
#'   (passed to [fp_data()]).
#' @param log_vars Character vector of variable names stored on the log scale
#'   (passed to [fp_data()]).
#' @param percs   Percentiles shown for continuous predictors
#'   (passed to [fp_data()], default `c(0.05, 0.25, 0.75, 0.95)`).
#' @param ci_level Confidence level (passed to [fp_data()], default `0.95`).
#' @param digits  Significant figures for formatted text
#'   (passed to [fp_data()], default `3`).
#' @param ci_width Spacer column width (passed to [fp_data()], default `34`).
#' @param xlim       Length-2 numeric; x-axis limits on the original
#'   (non-transformed) scale (default `c(0.1, 10)`).
#' @param vert_line  Numeric vector of additional vertical reference lines
#'   (default `c(0.8, 1.25)`).  Set to `NULL` to suppress.
#' @param ci_col     Colour string for CI points and bars (default dark red).
#' @param base_size  Base font size (default `10`).
#' @param x_trans    Axis transformation: `"log"` (default), `"log2"`,
#'   `"log10"`, or `"none"`.
#' @param ticks_at   Numeric vector of custom tick positions on the original
#'   x-axis scale.  `NULL` (default) uses automatic ticks.
#'
#' @return A `gtable` object of class `"forestplot"`.  Print it directly or
#'   draw it with `grid::grid.draw()` inside a graphics device.
#'
#' @seealso [fp_data()]
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   lung <- survival::lung
#'   lung$status <- lung$status - 1L
#'   lung$sex <- factor(lung$sex, levels = 1:2, labels = c("Male", "Female"))
#'   lung_cc <- lung[
#'     complete.cases(lung[, c("time", "status", "age", "sex")]), ]
#'
#'   fit <- survival::coxph(
#'     survival::Surv(time, status) ~ age + sex,
#'     data = lung_cc
#'   )
#'   fp_plot(fit, lung_cc,
#'           labels = c(age = "Age (years)", sex = "Sex"))
#' }
#'
#' @export
fp_plot <- function(
  model,
  data,
  labels    = NULL,
  log_vars  = NULL,
  percs     = c(0.05, 0.25, 0.75, 0.95),
  ci_level  = 0.95,
  digits    = 3L,
  ci_width  = 34L,
  xlim      = c(0.1, 10),
  vert_line = c(0.8, 1.25),
  ci_col    = "#C00000",
  base_size = 10,
  x_trans   = "log",
  ticks_at  = NULL
) {
  plot_data <- fp_data(
    model    = model,
    data     = data,
    labels   = labels,
    log_vars = log_vars,
    percs    = percs,
    ci_level = ci_level,
    digits   = digits,
    ci_width = ci_width
  )

  is_cox   <- "HR (95% CI)" %in% names(plot_data)
  est_col  <- if (is_cox) "HR (95% CI)" else "OR (95% CI)"
  arr_labs <- if (is_cox)
    c("Lower hazard", "Greater hazard")
  else
    c("Lower likelihood", "Greater likelihood")

  disp_data <- plot_data[, c("Predictor", est_col, " ", "P-value")]

  tm <- forestploter::forest_theme(
    base_size        = base_size,
    summary_col      = "black",
    arrow_label_just = "end",
    arrow_type       = "closed",
    ci_pch           = 16,
    ci_col           = ci_col,
    ci_fill          = ci_col,
    ci_alpha         = 1,
    ci_lty           = 1,
    ci_lwd           = 1.5,
    ci_Theight       = 0.2
  )

  fp <- forestploter::forest(
    disp_data,
    est       = plot_data$est,
    lower     = plot_data$lci,
    upper     = plot_data$uci,
    ci_column = 3L,
    ref_line  = 1,
    vert_line = vert_line,
    sizes     = 1,
    arrow_lab = arr_labs,
    x_trans   = x_trans,
    xlim      = xlim,
    ticks_at  = ticks_at,
    theme     = tm
  )

  hrows <- which(plot_data$is_header)
  if (length(hrows) > 0L) {
    fp <- forestploter::edit_plot(
      fp, col = seq_len(ncol(disp_data)), row = hrows,
      which = "background", gp = grid::gpar(fill = "ivory2")
    )
    fp <- forestploter::edit_plot(
      fp, col = 1L, row = hrows,
      gp = grid::gpar(fontface = "bold")
    )
  }

  sig_rows <- which(
    plot_data$is_header & !is.na(plot_data$pval) & plot_data$pval < 0.05
  )
  if (length(sig_rows) > 0L) {
    fp <- forestploter::edit_plot(fp, col = 4L, row = sig_rows,
                                  gp = grid::gpar(col = "red"))
  }

  forestploter::add_border(fp, part = "header", row = 1L,
                           gp = grid::gpar(lwd = 1))
}
