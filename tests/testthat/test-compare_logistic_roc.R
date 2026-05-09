# Tests for compare_logistic_roc()

make_roc_df <- function(n = 200, seed = 11L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 50, 10),
    bmi = rnorm(n, 27, 5),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  linpred <- -4 + 0.05 * df$age + 0.08 * df$bmi + 0.6 * (df$sex == "M")
  df$y <- rbinom(n, 1, plogis(linpred))
  df
}

test_that("compare_logistic_roc() returns a ggplot object", {
  df <- make_roc_df()
  fit_age <- glm(y ~ age, data = df, family = binomial)
  fit_full <- glm(y ~ age + bmi + sex, data = df, family = binomial)

  p <- compare_logistic_roc(list("Age only" = fit_age, "Full" = fit_full))

  expect_s3_class(p, "ggplot")
})

test_that("compare_logistic_roc() overlays one ROC curve per model", {
  df <- make_roc_df()
  fit_age <- glm(y ~ age, data = df, family = binomial)
  fit_full <- glm(y ~ age + bmi + sex, data = df, family = binomial)

  p <- compare_logistic_roc(list("Age only" = fit_age, "Full" = fit_full))
  auc_tbl <- plotROC::calc_auc(p)

  expect_equal(sort(as.character(auc_tbl$model)), c("Age only", "Full"))
  expect_equal(levels(p$data$model), c("Age only", "Full"))
})

test_that("compare_logistic_roc() can append AUC values to legend labels", {
  df <- make_roc_df()
  fit_age <- glm(y ~ age, data = df, family = binomial)
  fit_full <- glm(y ~ age + bmi + sex, data = df, family = binomial)

  p <- compare_logistic_roc(
    list("Age only" = fit_age, "Full" = fit_full),
    show_auc = TRUE,
    auc_digits = 2L
  )

  scale <- p$scales$get_scales("colour")
  labels <- scale$labels(c("Age only", "Full"))

  expect_match(labels[[1]], "Age only \\(AUC = ")
  expect_match(labels[[2]], "Full \\(AUC = ")
})

test_that("compare_logistic_roc() aligns models on shared observations", {
  df <- make_roc_df()
  df$bmi[1:20] <- NA_real_

  fit_age <- glm(y ~ age, data = df, family = binomial)
  fit_full <- glm(y ~ age + bmi + sex, data = df, family = binomial)

  p <- compare_logistic_roc(list("Age only" = fit_age, "Full" = fit_full))
  shared_n <- sum(stats::complete.cases(df[, c("y", "age", "bmi", "sex")]))

  expect_equal(nrow(p$data), 2L * shared_n)
})

test_that("compare_logistic_roc() rejects too few models", {
  df <- make_roc_df()
  fit_age <- glm(y ~ age, data = df, family = binomial)

  expect_error(
    compare_logistic_roc(list("Age only" = fit_age)),
    "at least two"
  )
})

test_that("compare_logistic_roc() rejects non-binomial fits", {
  df <- make_roc_df()
  fit_age <- glm(y ~ age, data = df, family = binomial)
  fit_pois <- glm(y ~ age, data = df, family = poisson)

  expect_error(
    compare_logistic_roc(list("Age only" = fit_age, "Poisson" = fit_pois)),
    "family = binomial"
  )
})

test_that("compare_logistic_roc() rejects mismatched outcomes", {
  df <- make_roc_df()
  df$y2 <- 1L - df$y
  fit_1 <- glm(y ~ age, data = df, family = binomial)
  fit_2 <- glm(y2 ~ age + bmi, data = df, family = binomial)

  expect_error(
    compare_logistic_roc(list("Model 1" = fit_1, "Model 2" = fit_2)),
    "same binary outcome"
  )
})
