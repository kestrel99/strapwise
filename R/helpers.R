# Internal helpers for building prediction grids and computing bootstrapped
# predicted probabilities.  These functions are not exported.

# .build_pred_frame ------------------------------------------------------------
#
# Build a prediction data frame for a single predictor, holding all other
# predictors at their median (numeric) or mode (factor / character).
#
# @param bl        A `boot_logistic` object.
# @param term_name Character scalar: base variable name as it appears in the
#   original data (e.g. `"age"`, not `"log(age)"`).
# @param n_grid    Number of equally-spaced grid points for continuous
#   predictors.  Ignored for factors.
#
# @return A data frame suitable for passing to `model.matrix()`.
.build_pred_frame <- function(bl, term_name, n_grid = 200L) {
  data <- bl$data
  fit <- bl$fit

  # Derive response and predictor names from the formula so the lookup is
  # robust regardless of column order and inline transformations.
  resp_name <- as.character(formula(fit)[[2]])
  pred_cols <- intersect(all.vars(formula(fit)[[3]]), names(data))

  focal_col <- data[[term_name]]
  is_factor <- is.factor(focal_col) || is.character(focal_col)

  # Reference values for all non-focal predictors
  ref_vals <- lapply(stats::setNames(pred_cols, pred_cols), function(v) {
    col <- data[[v]]
    if (v == term_name) {
      return(NULL)
    }
    if (is.numeric(col)) {
      return(stats::median(col, na.rm = TRUE))
    }
    tbl <- sort(table(col), decreasing = TRUE)
    factor(names(tbl)[1L], levels = levels(factor(col)))
  })

  # Focal grid
  focal_grid <- if (is_factor) {
    factor(levels(factor(focal_col)), levels = levels(factor(focal_col)))
  } else {
    seq(
      min(focal_col, na.rm = TRUE),
      max(focal_col, na.rm = TRUE),
      length.out = n_grid
    )
  }

  # Assemble prediction frame
  non_focal <- pred_cols[pred_cols != term_name]
  pred_list <- c(
    stats::setNames(list(focal_grid), term_name),
    lapply(
      stats::setNames(non_focal, non_focal),
      function(v) rep(ref_vals[[v]], length(focal_grid))
    )
  )
  pred_df <- as.data.frame(pred_list, stringsAsFactors = FALSE)

  # Restore factor levels from original data
  for (v in pred_cols) {
    if (is.factor(data[[v]])) {
      pred_df[[v]] <- factor(pred_df[[v]], levels = levels(data[[v]]))
    }
  }

  # Append a dummy response column (required by model.matrix for some formulae)
  pred_df[[resp_name]] <- data[[resp_name]][1L]

  pred_df
}


# .compute_preds ---------------------------------------------------------------
#
# Compute point-estimate predicted probabilities and percentile bootstrap CIs
# for a prediction data frame produced by `.build_pred_frame()`.
#
# @param bl      A `boot_logistic` object.
# @param pred_df A data frame from `.build_pred_frame()`.
#
# @return A data frame with columns `prob`, `ci_lo`, and `ci_hi`.
.compute_preds <- function(bl, pred_df) {
  fit <- bl$fit
  alpha <- 1 - bl$conf_level
  link_inv <- fit$family$linkinv

  X_pred <- stats::model.matrix(stats::terms(fit), data = pred_df)

  point_est <- as.numeric(link_inv(X_pred %*% coef(fit)))
  boot_preds <- apply(
    bl$boot_coefs,
    1,
    function(b) as.numeric(link_inv(X_pred %*% b))
  )

  ci_lo <- apply(
    boot_preds,
    1,
    stats::quantile,
    probs = alpha / 2,
    na.rm = TRUE
  )
  ci_hi <- apply(
    boot_preds,
    1,
    stats::quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )

  data.frame(prob = point_est, ci_lo = ci_lo, ci_hi = ci_hi)
}


# .obs_proportions -------------------------------------------------------------
#
# Compute observed response proportions (and Wilson score CIs) for a continuous
# predictor binned into quantile strata, or for each level of a categorical
# predictor.  Optionally stratified by a grouping variable.
#
# For continuous `x` the observations are split into quantile bins.  When
# `group` is NULL the requested `n_groups` is applied to the full data.  When
# `group` is supplied, binning is performed *per group* using each group's own
# quantile breaks, so that strata reflect the within-group distribution rather
# than the pooled one.  If a group does not have enough distinct values to
# support `n_groups` bins, the bin count is reduced by one and retried, down to
# a minimum of 2 bins.  A group is omitted only if it has fewer than 2
# observations (so no split is possible at all).
#
# For factor/character `x` each level is its own stratum; binning and
# `n_groups` are ignored.
#
# The Wilson score interval is used for CIs on proportions because it maintains
# near-nominal coverage even for small `n` or extreme proportions, unlike the
# normal approximation.
#
# @param x          Numeric or factor/character predictor vector.
# @param y          Numeric 0/1 response vector (same length as `x`).
# @param n_groups   Integer >= 2; maximum number of quantile strata for
#   continuous `x`.  The actual bin count may be lower for a given group if
#   there are insufficient distinct values.  Ignored for factor/character `x`.
# @param conf_level Confidence level for Wilson CIs.
# @param group      Optional factor/character grouping vector (same length as
#   `x`).  When not `NULL`, proportions are computed per (group x stratum) with
#   per-group quantile breaks.  A `grp` column is added to the result.
#
# @return A data frame with columns `x`, `obs`, `ci_lo`, `ci_hi`, `n`, and
#   (when `group` is supplied) `grp`.  Returns `NULL` only when there are
#   fewer than 2 complete observations in total.
.obs_proportions <- function(
  x,
  y,
  n_groups = 4L,
  conf_level = 0.95,
  group = NULL
) {
  # ---- Complete-case filtering -----------------------------------------------
  ok_vars <- if (is.null(group)) list(x, y) else list(x, y, group)
  ok <- Reduce("&", lapply(ok_vars, function(v) !is.na(v)))
  x <- x[ok]
  y <- as.numeric(y[ok])
  if (!is.null(group)) {
    group <- as.character(group[ok])
  }

  if (length(x) < 2L) {
    return(NULL)
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  is_factor <- is.factor(x) || is.character(x)

  # ---- Wilson CI for a single proportion -------------------------------------
  .wilson <- function(p, n) {
    denom <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denom
    margin <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
    c(lo = max(0, center - margin), hi = min(1, center + margin))
  }

  # ---- Derive breaks for xs, trying n_req bins and falling back downward -----
  # Returns a valid breaks vector (>= 3 unique values), or NULL when all values
  # in xs are identical (meaning no split is possible at any bin count).
  .make_breaks <- function(xs, n_req) {
    while (n_req >= 2L) {
      probs <- seq(0, 1, length.out = n_req + 1L)
      breaks <- unique(stats::quantile(xs, probs = probs, na.rm = TRUE))
      if (length(breaks) >= 3L) {
        return(breaks)
      }
      n_req <- n_req - 1L
    }
    NULL # all values identical  -- 1-bin fallback handled in .summarise_one
  }

  # ---- Summarise one (xs, ys) subset into stratum rows -----------------------
  .summarise_one <- function(xs, ys) {
    if (length(xs) < 1L) {
      return(NULL)
    }

    if (is_factor) {
      lvls <- levels(factor(x)) # global levels for consistent ordering
      strata <- factor(xs, levels = lvls)
      strata_ids <- lvls
      get_pos <- function(g, xsub) g
    } else {
      breaks <- .make_breaks(xs, n_groups)
      if (is.null(breaks)) {
        # All values identical: treat the whole subset as a single stratum.
        n <- length(ys)
        p <- mean(ys, na.rm = TRUE)
        ci <- .wilson(p, n)
        return(data.frame(
          x = mean(xs, na.rm = TRUE),
          obs = p,
          ci_lo = ci[["lo"]],
          ci_hi = ci[["hi"]],
          n = n,
          stringsAsFactors = FALSE
        ))
      }
      strata <- cut(xs, breaks = breaks, include.lowest = TRUE, labels = FALSE)
      strata_ids <- sort(unique(strata[!is.na(strata)]))
      get_pos <- function(g, xsub) mean(xsub, na.rm = TRUE)
    }

    rows <- lapply(strata_ids, function(g) {
      idx <- if (is_factor) strata == g else (!is.na(strata) & strata == g)
      if (!any(idx)) {
        return(NULL)
      }
      xg <- xs[idx]
      yg <- ys[idx]
      n <- length(yg)
      if (n < 1L) {
        return(NULL)
      }
      p <- mean(yg, na.rm = TRUE)
      ci <- .wilson(p, n)
      data.frame(
        x = get_pos(g, xg),
        obs = p,
        ci_lo = ci[["lo"]],
        ci_hi = ci[["hi"]],
        n = n,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }

  # ---- Dispatch: ungrouped or grouped ----------------------------------------
  if (is.null(group)) {
    .summarise_one(x, y)
  } else {
    grp_lvls <- unique(group)
    parts <- lapply(grp_lvls, function(g) {
      idx <- group == g
      res <- .summarise_one(x[idx], y[idx])
      if (is.null(res) || nrow(res) == 0L) {
        return(NULL)
      }
      res$grp <- g
      res
    })
    out <- do.call(rbind, Filter(Negate(is.null), parts))
    if (is.null(out) || nrow(out) == 0L) {
      return(NULL)
    }
    out
  }
}
