# Tests for plot_logistic_curves().

# ---- Input validation -------------------------------------------------------

test_that("plot_logistic_curves() rejects non-boot_logistic input", {
  expect_error(plot_logistic_curves(FIT), "boot_logistic")
  expect_error(plot_logistic_curves(list()), "boot_logistic")
})

# ---- Return type and structure ----------------------------------------------

test_that("plot_logistic_curves() returns a named list by default", {
  plots <- plot_logistic_curves(BL)
  expect_type(plots, "list")
})

test_that("list has one entry per predictor", {
  plots <- plot_logistic_curves(BL)
  pred_terms <- intersect(all.vars(formula(FIT)[[3]]), names(DF))
  expect_equal(length(plots), length(pred_terms))
})

test_that("list names are base variable names (not inline expressions)", {
  plots <- plot_logistic_curves(BL)
  pred_terms <- intersect(all.vars(formula(FIT)[[3]]), names(DF))
  expect_equal(sort(names(plots)), sort(pred_terms))
})

test_that("each element is a ggplot object", {
  plots <- plot_logistic_curves(BL)
  for (nm in names(plots)) {
    expect_s3_class(plots[[nm]], "ggplot")
  }
})

# ---- Individual predictor plots ---------------------------------------------

test_that("continuous predictor plot has geom_line and geom_ribbon", {
  p <- plot_logistic_curves(BL)[["age"]]
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomLine" %in% geoms)
  expect_true("GeomRibbon" %in% geoms)
})

test_that("continuous predictor plot has raw jitter layer by default", {
  p <- plot_logistic_curves(BL)[["age"]]
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomPoint" %in% geoms)
})

test_that("continuous predictor plot has no point layer when raw_data = FALSE", {
  plots <- plot_logistic_curves(BL, raw_data = FALSE, obs_groups = NULL)
  p <- plots[["age"]]
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_false("GeomPoint" %in% geoms)
})

test_that("categorical predictor plot has error bar and point layers", {
  p <- plot_logistic_curves(BL)[["trt"]]
  geoms <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomErrorbar" %in% geoms)
  expect_true("GeomPoint" %in% geoms)
})

# ---- n_grid controls curve resolution  -------------------------------------

test_that("n_grid controls number of fitted curve points", {
  p <- plot_logistic_curves(BL, n_grid = 25L)[["age"]]
  bd <- ggplot2::ggplot_build(p)
  # First built-data slot corresponds to the ribbon (first layer); should
  # have exactly n_grid = 25 rows.
  expect_equal(nrow(bd$data[[1]]), 25L)
})

# ---- Inline transformation: plot accessible by base variable name -----------

test_that("log() transformation: plot keyed by base name 'age'", {
  fit_log <- glm(y ~ log(age) + sex, data = DF, family = binomial)
  bl_log <- bootstrap_logistic(fit_log, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_log)
  expect_true("age" %in% names(plots))
  expect_false("log(age)" %in% names(plots))
  expect_s3_class(plots[["age"]], "ggplot")
})

# ---- Aesthetic parameters ---------------------------------------------------

test_that("custom palette is accepted without error", {
  expect_no_error(
    plot_logistic_curves(BL, palette = c("steelblue", "tomato", "forestgreen"))
  )
})

test_that("custom base_theme is applied without error", {
  expect_no_error(
    plot_logistic_curves(
      BL,
      base_theme = ggplot2::theme_minimal(base_size = 14)
    )
  )
})

test_that("free_y = TRUE removes fixed y-axis limits", {
  p_fixed <- plot_logistic_curves(BL, free_y = FALSE)[["age"]]
  p_free <- plot_logistic_curves(BL, free_y = TRUE)[["age"]]
  lims_fixed <- ggplot2::layer_scales(p_fixed)$y$limits
  lims_free <- ggplot2::layer_scales(p_free)$y$limits
  expect_equal(lims_fixed, c(0, 1))
  expect_null(lims_free)
})

test_that("raw_alpha, raw_point_size, jitter_height arguments accepted", {
  expect_no_error(
    plot_logistic_curves(
      BL,
      raw_alpha = 0.5,
      raw_point_size = 2,
      jitter_height = 0.05
    )
  )
})

# ---- combine argument -------------------------------------------------------

test_that("combine = FALSE returns a list (default)", {
  result <- plot_logistic_curves(BL, combine = FALSE)
  expect_type(result, "list")
})

test_that("combine = TRUE returns patchwork or list with message", {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    result <- plot_logistic_curves(BL, combine = TRUE)
    expect_true(inherits(result, "patchwork") || is.list(result))
  } else {
    expect_message(plot_logistic_curves(BL, combine = TRUE), "patchwork")
  }
})

# ---- Alternative link functions ---------------------------------------------

test_that("plots render for probit link", {
  fit_p <- glm(y ~ age + sex, data = DF, family = binomial("probit"))
  bl_p <- bootstrap_logistic(fit_p, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_p)
  expect_s3_class(plots[["age"]], "ggplot")
  expect_s3_class(plots[["sex"]], "ggplot")
})

test_that("plots render for cloglog link", {
  fit_c <- glm(y ~ age + sex, data = DF, family = binomial("cloglog"))
  bl_c <- bootstrap_logistic(fit_c, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_c)
  expect_s3_class(plots[["age"]], "ggplot")
})

# ---- Intercept-only and single-predictor models -----------------------------

test_that("intercept-only model returns empty list without error", {
  fit_null <- glm(y ~ 1, data = DF, family = binomial)
  bl_null <- bootstrap_logistic(fit_null, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_null)
  expect_equal(length(plots), 0L)
})

test_that("single continuous predictor produces one plot", {
  fit_1 <- glm(y ~ age, data = DF, family = binomial)
  bl_1 <- bootstrap_logistic(fit_1, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_1)
  expect_equal(length(plots), 1L)
  expect_s3_class(plots[["age"]], "ggplot")
})

test_that("single categorical predictor produces one plot", {
  fit_cat <- glm(y ~ trt, data = DF, family = binomial)
  bl_cat <- bootstrap_logistic(fit_cat, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_cat)
  expect_equal(length(plots), 1L)
  expect_s3_class(plots[["trt"]], "ggplot")
})

# ---- ggplot build (catches runtime rendering errors) -----------------------

test_that("all plots build without error", {
  plots <- plot_logistic_curves(BL)
  for (nm in names(plots)) {
    expect_no_error(ggplot2::ggplot_build(plots[[nm]]))
  }
})

test_that("plots build without error when raw_data = FALSE", {
  plots <- plot_logistic_curves(BL, raw_data = FALSE)
  for (nm in names(plots)) {
    expect_no_error(ggplot2::ggplot_build(plots[[nm]]))
  }
})

test_that("log-transformed plot builds without error", {
  fit_log <- glm(y ~ log(age) + sex, data = DF, family = binomial)
  bl_log <- bootstrap_logistic(fit_log, n_boot = 100L, seed = 1L)
  plots <- plot_logistic_curves(bl_log)
  for (nm in names(plots)) {
    expect_no_error(ggplot2::ggplot_build(plots[[nm]]))
  }
})
