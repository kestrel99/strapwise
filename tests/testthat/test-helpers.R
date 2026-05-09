# Tests for internal helpers: .build_pred_frame() and .compute_preds()

# ---- .build_pred_frame() ----------------------------------------------------

test_that(".build_pred_frame() returns a data frame", {
  pf <- strapwise:::.build_pred_frame(BL, "age")
  expect_s3_class(pf, "data.frame")
})

test_that(".build_pred_frame() grid length matches n_grid for continuous", {
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = 50L)
  expect_equal(nrow(pf), 50L)
})

test_that(".build_pred_frame() has one row per level for categorical", {
  pf <- strapwise:::.build_pred_frame(BL, "trt", n_grid = 1L)
  n_levels <- nlevels(DF$trt)
  expect_equal(nrow(pf), n_levels)
})

test_that(".build_pred_frame() focal column spans observed range", {
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = 100L)
  expect_gte(min(pf$age), min(DF$age) - .Machine$double.eps^0.5)
  expect_lte(max(pf$age), max(DF$age) + .Machine$double.eps^0.5)
})

test_that(".build_pred_frame() fixes non-focal numeric at median", {
  med_bmi <- median(DF$bmi)
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = 20L)
  expect_true(all(abs(pf$bmi - med_bmi) < 1e-10))
})

test_that(".build_pred_frame() fixes non-focal factor at mode", {
  mode_sex <- names(sort(table(DF$sex), decreasing = TRUE))[1]
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = 20L)
  expect_true(all(as.character(pf$sex) == mode_sex))
})

test_that(".build_pred_frame() preserves factor levels", {
  pf <- strapwise:::.build_pred_frame(BL, "trt", n_grid = 1L)
  expect_equal(levels(pf$trt), levels(DF$trt))
})

test_that(".build_pred_frame() includes the response column", {
  resp <- as.character(formula(FIT)[[2]])
  pf <- strapwise:::.build_pred_frame(BL, "age")
  expect_true(resp %in% names(pf))
})

test_that(".build_pred_frame() works for each predictor in the model", {
  for (term in c("age", "bmi", "sex", "trt")) {
    pf <- strapwise:::.build_pred_frame(BL, term, n_grid = 20L)
    expect_s3_class(pf, "data.frame")
    expect_true(term %in% names(pf))
  }
})

test_that(".build_pred_frame() works with inline transformation (log)", {
  fit_log <- glm(y ~ log(age) + sex, data = DF, family = binomial)
  bl_log <- bootstrap_logistic(fit_log, n_boot = 100L, seed = 1L)
  pf <- strapwise:::.build_pred_frame(bl_log, "age", n_grid = 30L)
  expect_equal(nrow(pf), 30L)
  expect_true("age" %in% names(pf))
})

# ---- .compute_preds() -------------------------------------------------------

test_that(".compute_preds() returns a data frame with correct columns", {
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = 50L)
  res <- strapwise:::.compute_preds(BL, pf)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("prob", "ci_lo", "ci_hi"))
})

test_that(".compute_preds() returns one row per prediction grid point", {
  n_grid <- 40L
  pf <- strapwise:::.build_pred_frame(BL, "age", n_grid = n_grid)
  res <- strapwise:::.compute_preds(BL, pf)
  expect_equal(nrow(res), n_grid)
})

test_that("predicted probabilities are in [0, 1]", {
  pf <- strapwise:::.build_pred_frame(BL, "age")
  res <- strapwise:::.compute_preds(BL, pf)
  expect_true(all(res$prob >= 0 & res$prob <= 1))
  expect_true(all(res$ci_lo >= 0 & res$ci_lo <= 1))
  expect_true(all(res$ci_hi >= 0 & res$ci_hi <= 1))
})

test_that("ci_lo <= prob <= ci_hi for all grid points", {
  pf <- strapwise:::.build_pred_frame(BL, "age")
  res <- strapwise:::.compute_preds(BL, pf)
  expect_true(all(res$ci_lo <= res$prob + 1e-10))
  expect_true(all(res$prob <= res$ci_hi + 1e-10))
})

test_that("CI width is positive for all grid points", {
  pf <- strapwise:::.build_pred_frame(BL, "age")
  res <- strapwise:::.compute_preds(BL, pf)
  expect_true(all((res$ci_hi - res$ci_lo) > 0))
})

test_that(".compute_preds() respects linkinv for probit link", {
  fit_p <- glm(y ~ age + sex, data = DF, family = binomial("probit"))
  bl_p <- bootstrap_logistic(fit_p, n_boot = 100L, seed = 1L)
  pf <- strapwise:::.build_pred_frame(bl_p, "age", n_grid = 30L)
  res <- strapwise:::.compute_preds(bl_p, pf)
  expect_true(all(res$prob >= 0 & res$prob <= 1))
})

test_that(".compute_preds() works for categorical predictor grid", {
  pf <- strapwise:::.build_pred_frame(BL, "trt", n_grid = 1L)
  res <- strapwise:::.compute_preds(BL, pf)
  expect_equal(nrow(res), nlevels(DF$trt))
  expect_true(all(res$prob >= 0 & res$prob <= 1))
})

test_that(".compute_preds() 90% CI is narrower than 95% CI", {
  bl90 <- bootstrap_logistic(FIT, n_boot = 150L, conf_level = 0.90, seed = 1L)
  bl95 <- bootstrap_logistic(FIT, n_boot = 150L, conf_level = 0.95, seed = 1L)
  pf <- strapwise:::.build_pred_frame(bl90, "age", n_grid = 20L)
  r90 <- strapwise:::.compute_preds(bl90, pf)
  r95 <- strapwise:::.compute_preds(bl95, pf)
  width90 <- mean(r90$ci_hi - r90$ci_lo)
  width95 <- mean(r95$ci_hi - r95$ci_lo)
  expect_lt(width90, width95)
})
