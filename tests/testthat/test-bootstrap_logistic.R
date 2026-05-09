# Tests for bootstrap_logistic() and its S3 methods.

# ---- Input validation -------------------------------------------------------

test_that("bootstrap_logistic() rejects non-glm input", {
  expect_error(
    bootstrap_logistic(lm(y ~ age, data = DF)),
    "`fit` must be a fitted glm object"
  )
})

test_that("bootstrap_logistic() rejects non-binomial family", {
  fit_pois <- glm(y ~ age, data = DF, family = poisson)
  expect_error(bootstrap_logistic(fit_pois), "family = binomial")
})

test_that("bootstrap_logistic() rejects out-of-range conf_level", {
  expect_error(bootstrap_logistic(FIT, conf_level = 0), "conf_level")
  expect_error(bootstrap_logistic(FIT, conf_level = 1), "conf_level")
  expect_error(bootstrap_logistic(FIT, conf_level = 1.5), "conf_level")
  expect_error(bootstrap_logistic(FIT, conf_level = -0.5), "conf_level")
})

# ---- Return structure -------------------------------------------------------

test_that("bootstrap_logistic() returns a boot_logistic object", {
  expect_s3_class(BL, "boot_logistic")
})

test_that("boot_logistic object has all required fields", {
  expect_named(
    BL,
    c("fit", "boot_coefs", "data", "formula", "conf_level", "n_boot"),
    ignore.order = TRUE
  )
})

test_that("boot_coefs has correct number of columns", {
  expect_equal(ncol(BL$boot_coefs), length(coef(FIT)))
})

test_that("boot_coefs column names match model coefficient names", {
  expect_equal(colnames(BL$boot_coefs), names(coef(FIT)))
})

test_that("n_boot equals actual number of rows in boot_coefs", {
  expect_equal(BL$n_boot, nrow(BL$boot_coefs))
})

test_that("n_boot is at most the requested number of replicates", {
  expect_lte(BL$n_boot, 150L)
})

test_that("conf_level is stored correctly", {
  bl90 <- bootstrap_logistic(FIT, n_boot = 100L, conf_level = 0.90, seed = 1L)
  expect_equal(bl90$conf_level, 0.90)
})

test_that("formula field matches the original model formula", {
  expect_equal(deparse(BL$formula), deparse(formula(FIT)))
})

test_that("data field is a data frame with the expected columns", {
  expect_s3_class(BL$data, "data.frame")
  expect_true(all(c("age", "bmi", "sex", "trt", "y") %in% names(BL$data)))
})

test_that("fit field is the original glm object", {
  expect_identical(BL$fit, FIT)
})

# ---- boot_coefs integrity ---------------------------------------------------

test_that("boot_coefs contains no NA after removal step", {
  expect_false(anyNA(BL$boot_coefs))
})

test_that("bootstrapped intercepts are finite", {
  expect_true(all(is.finite(BL$boot_coefs[, "(Intercept)"])))
})

test_that("bootstrap SE is positive for all coefficients", {
  sds <- apply(BL$boot_coefs, 2, sd)
  expect_true(all(sds > 0))
})

# ---- Reproducibility --------------------------------------------------------

test_that("same seed produces identical bootstrap draws", {
  bl_a <- bootstrap_logistic(FIT, n_boot = 100L, seed = 99L)
  bl_b <- bootstrap_logistic(FIT, n_boot = 100L, seed = 99L)
  expect_equal(bl_a$boot_coefs, bl_b$boot_coefs)
})

test_that("different seeds produce different bootstrap draws", {
  bl_a <- bootstrap_logistic(FIT, n_boot = 100L, seed = 1L)
  bl_b <- bootstrap_logistic(FIT, n_boot = 100L, seed = 2L)
  expect_false(identical(bl_a$boot_coefs, bl_b$boot_coefs))
})

# ---- Alternative link functions ---------------------------------------------

test_that("bootstrap_logistic() works with probit link", {
  fit_p <- glm(y ~ age + sex, data = DF, family = binomial("probit"))
  bl_p <- bootstrap_logistic(fit_p, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_p, "boot_logistic")
  expect_equal(bl_p$fit$family$link, "probit")
})

test_that("bootstrap_logistic() works with cloglog link", {
  fit_c <- glm(y ~ age + sex, data = DF, family = binomial("cloglog"))
  bl_c <- bootstrap_logistic(fit_c, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_c, "boot_logistic")
})

test_that("bootstrap_logistic() works with cauchit link", {
  fit_ca <- glm(y ~ age, data = DF, family = binomial("cauchit"))
  bl_ca <- bootstrap_logistic(fit_ca, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_ca, "boot_logistic")
})

# ---- Inline formula transformations ----------------------------------------

test_that("log() transformation: boot_coefs names and data columns correct", {
  fit_log <- glm(y ~ log(age) + sex, data = DF, family = binomial)
  bl_log <- bootstrap_logistic(fit_log, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_log, "boot_logistic")
  expect_true("log(age)" %in% colnames(bl_log$boot_coefs))
  expect_true("age" %in% names(bl_log$data))
  expect_false("log(age)" %in% names(bl_log$data))
})

test_that("sqrt() transformation works", {
  df2 <- DF
  df2$bmi <- abs(df2$bmi) + 1
  fit_sqrt <- glm(y ~ sqrt(bmi) + sex, data = df2, family = binomial)
  bl_sqrt <- bootstrap_logistic(fit_sqrt, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_sqrt, "boot_logistic")
  expect_true("bmi" %in% names(bl_sqrt$data))
})

test_that("I() wrapper works", {
  fit_i <- glm(y ~ I(age^2) + sex, data = DF, family = binomial)
  bl_i <- bootstrap_logistic(fit_i, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_i, "boot_logistic")
})

# ---- Edge-case models -------------------------------------------------------

test_that("intercept-only model works", {
  fit_null <- glm(y ~ 1, data = DF, family = binomial)
  bl_null <- bootstrap_logistic(fit_null, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_null, "boot_logistic")
  expect_equal(ncol(bl_null$boot_coefs), 1L)
})

test_that("single continuous predictor works", {
  fit_1 <- glm(y ~ age, data = DF, family = binomial)
  bl_1 <- bootstrap_logistic(fit_1, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_1, "boot_logistic")
})

test_that("single categorical predictor works", {
  fit_cat <- glm(y ~ trt, data = DF, family = binomial)
  bl_cat <- bootstrap_logistic(fit_cat, n_boot = 100L, seed = 1L)
  expect_s3_class(bl_cat, "boot_logistic")
  expect_equal(ncol(bl_cat$boot_coefs), 3L) # intercept + 2 contrasts
})

# ---- Warning when too few replicates succeed --------------------------------

test_that("warning fires when successful replicates drop below 100", {
  # Near-perfect separation on a tiny dataset forces most resamples to have
  # extreme / NA coefficients.
  set.seed(42)
  df_sep <- data.frame(x = c(1:5, 20:24), y = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))
  fit_sep <- suppressWarnings(glm(y ~ x, data = df_sep, family = binomial))
  expect_warning(
    suppressMessages(bootstrap_logistic(fit_sep, n_boot = 50L, seed = 1L)),
    regexp = "replicates succeeded"
  )
})

test_that("warning-generating bootstrap refits are skipped and counted", {
  set.seed(42)
  df_sep <- data.frame(x = c(1:5, 20:24), y = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))
  fit_sep <- suppressWarnings(glm(y ~ x, data = df_sep, family = binomial))

  expect_message(
    suppressWarnings(bootstrap_logistic(fit_sep, n_boot = 20L, seed = 1L)),
    regexp = "skipped due to errors or warnings"
  )
})

# ---- S3 print method --------------------------------------------------------

test_that("print.boot_logistic() runs without error", {
  expect_output(print(BL))
})

test_that("print.boot_logistic() returns the object invisibly", {
  out <- withVisible(print(BL))
  expect_false(out$visible)
  expect_identical(out$value, BL)
})

test_that("print output contains key fields", {
  out <- capture.output(print(BL))
  expect_true(any(grepl("Formula", out)))
  expect_true(any(grepl("Link", out)))
  expect_true(any(grepl("Replicates", out)))
  expect_true(any(grepl("CI level", out)))
})

# ---- S3 summary method ------------------------------------------------------

test_that("summary.boot_logistic() returns a data frame invisibly", {
  result <- withVisible(summary(BL))
  expect_false(result$visible)
  expect_s3_class(result$value, "data.frame")
})

test_that("summary table has correct dimensions", {
  tbl <- suppressMessages(summary(BL))
  expect_equal(nrow(tbl), length(coef(FIT)))
  expect_equal(ncol(tbl), 5L)
})

test_that("summary column names reflect conf_level", {
  tbl <- suppressMessages(summary(BL))
  ci_pct <- round(BL$conf_level * 100)
  expect_true(any(grepl(as.character(ci_pct), names(tbl))))
})

test_that("summary CI_lo <= Estimate <= CI_hi for all coefficients", {
  tbl <- suppressMessages(summary(BL))
  expect_true(all(tbl[[4]] <= tbl$Estimate))
  expect_true(all(tbl$Estimate <= tbl[[5]]))
})

test_that("summary row names match coefficient names", {
  tbl <- suppressMessages(summary(BL))
  expect_equal(rownames(tbl), names(coef(FIT)))
})

test_that("summary Boot_SE values are non-negative", {
  tbl <- suppressMessages(summary(BL))
  expect_true(all(tbl$Boot_SE >= 0))
})

test_that("summary conf_level = 0.90 produces correct column label", {
  bl90 <- bootstrap_logistic(FIT, n_boot = 100L, conf_level = 0.90, seed = 1L)
  tbl <- suppressMessages(summary(bl90))
  expect_true(any(grepl("90", names(tbl))))
})
