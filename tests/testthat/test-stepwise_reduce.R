# Tests for stepwise_reduce() and its S3 methods.
#
# Correctness tests compare against drop1() at alpha = 0.05 to approximate
# a realistic reduction scenario.  drop1() is the reference implementation;
# manual LRT mirrors the internals of stepwise_reduce().

library(survival)

# ---- Local fixtures ----------------------------------------------------------

# Logistic: age is significant at 0.05; bmi, sbp, crp are not.
make_log_df <- function(seed = 42L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 50, 12),
    bmi = rnorm(n, 27,  5),
    sbp = rnorm(n, 130, 20),
    crp = rexp(n, 0.5)
  )
  df$y <- rbinom(n, 1, plogis(-2 + 0.04 * df$age))
  df
}

# Logistic with 4-level factor: age significant, group + score not.
make_fac_df <- function(seed = 7L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age   = rnorm(n, 50, 12),
    group = factor(sample(c("A","B","C","D"), n, replace = TRUE)),
    score = rnorm(n)
  )
  df$y <- rbinom(n, 1, plogis(-1 + 0.05 * df$age))
  df
}

# Logistic interaction model.
make_int_df <- function(seed = 13L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 50, 12),
    sex = factor(sample(c("F","M"), n, replace = TRUE))
  )
  df$y <- rbinom(n, 1, plogis(-2 + 0.04 * df$age))
  df
}

# Cox: age mildly prognostic; bmi, score, noise not.
make_cox_df <- function(seed = 99L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age   = rnorm(n, 55, 12),
    bmi   = rnorm(n, 27,  5),
    score = rnorm(n),
    noise = rnorm(n)
  )
  df$time   <- rexp(n, exp(0.04 * df$age / 10))
  df$status <- rbinom(n, 1, 0.75)
  df
}

# Cox with 3-level factor.
make_cox_fac_df <- function(seed = 55L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age   = rnorm(n, 55, 12),
    ecog  = factor(sample(0:2, n, replace = TRUE)),
    noise = rnorm(n)
  )
  df$time   <- rexp(n, exp(0.03 * df$age / 10))
  df$status <- rbinom(n, 1, 0.70)
  df
}

LOG_DF  <- make_log_df()
FAC_DF  <- make_fac_df()
INT_DF  <- make_int_df()
COX_DF  <- make_cox_df()
COXF_DF <- make_cox_fac_df()

LOG_FIT  <- glm(y ~ age + bmi + sbp + crp,       data = LOG_DF,  family = binomial)
FAC_FIT  <- glm(y ~ age + group + score,           data = FAC_DF,  family = binomial)
INT_FIT  <- glm(y ~ age * sex,                     data = INT_DF,  family = binomial)
COX_FIT  <- coxph(Surv(time, status) ~ age + bmi + score + noise, data = COX_DF)
COXF_FIT <- coxph(Surv(time, status) ~ age + ecog + noise,        data = COXF_DF)

# Reductions (verbose = FALSE throughout)
LOG_SR  <- stepwise_reduce(LOG_FIT,  alpha = 0.05, verbose = FALSE)
FAC_SR  <- stepwise_reduce(FAC_FIT,  alpha = 0.05, verbose = FALSE)
INT_SR  <- stepwise_reduce(INT_FIT,  alpha = 0.05, verbose = FALSE)
COX_SR  <- stepwise_reduce(COX_FIT,  alpha = 0.05, verbose = FALSE)
COXF_SR <- stepwise_reduce(COXF_FIT, alpha = 0.05, verbose = FALSE)

# ---- Helpers for drop1-based comparison -------------------------------------

# Extract p-values from drop1() for all hierarchically testable terms.
drop1_pvals <- function(m) {
  is_cox <- inherits(m, "coxph")
  d1  <- if (is_cox) drop1(m, test = "Chisq") else drop1(m, test = "LRT")
  pv  <- d1[rownames(d1) != "<none>", , drop = FALSE]
  col <- grep("^Pr\\(", colnames(pv), value = TRUE)
  v   <- pv[[col[1L]]]
  names(v) <- rownames(pv)
  v
}

# Manual LRT p-values — mirrors stepwise_reduce internals.
manual_lrt_pvals <- function(m, dat) {
  is_cox <- inherits(m, "coxph")
  ties   <- if (is_cox) m$method else NULL
  trms   <- attr(stats::terms(m), "term.labels")
  vapply(setNames(trms, trms), function(trm) {
    new_f <- update.formula(formula(m), paste(". ~ . -", trm))
    red <- if (is_cox) {
      tryCatch(coxph(new_f, data = dat, ties = ties), error = function(e) NULL)
    } else {
      tryCatch(glm(new_f, data = dat, family = m$family), error = function(e) NULL)
    }
    if (is.null(red)) return(NA_real_)
    df_d <- max(1L, length(coef(m)) - length(coef(red)))
    lrt  <- max(0, -2 * (as.numeric(logLik(red)) - as.numeric(logLik(m))))
    pchisq(lrt, df = df_d, lower.tail = FALSE)
  }, numeric(1L))
}

# Term that drop1 would drop first (highest p-value above alpha = 0.05).
drop1_first_drop <- function(m) {
  pv  <- drop1_pvals(m)
  drp <- pv[pv >= 0.05]
  if (length(drp) == 0L) return(NA_character_)
  names(which.max(drp))
}

# ===========================================================================
# INPUT VALIDATION
# ===========================================================================

test_that("stepwise_reduce() rejects non-glm / non-coxph input", {
  expect_error(stepwise_reduce(lm(y ~ age, data = LOG_DF)),
               "`fit` must be a binomial glm or a coxph object")
})

test_that("stepwise_reduce() rejects non-binomial glm", {
  fit_p <- glm(y ~ age, data = LOG_DF, family = poisson)
  expect_error(stepwise_reduce(fit_p),
               "`fit` must be a binomial glm or a coxph object")
})

test_that("stepwise_reduce() rejects alpha outside (0, 1)", {
  expect_error(stepwise_reduce(LOG_FIT, alpha = 0),   "`alpha`")
  expect_error(stepwise_reduce(LOG_FIT, alpha = 1),   "`alpha`")
  expect_error(stepwise_reduce(LOG_FIT, alpha = -0.1), "`alpha`")
  expect_error(stepwise_reduce(LOG_FIT, alpha = 1.5),  "`alpha`")
})

test_that("stepwise_reduce() rejects non-character keep_vars", {
  expect_error(stepwise_reduce(LOG_FIT, keep_vars = 1L), "`keep_vars`")
})

test_that("stepwise_reduce() warns on keep_vars not in the model", {
  expect_warning(
    stepwise_reduce(LOG_FIT, keep_vars = "zzz", verbose = FALSE),
    "not in the model"
  )
})

# ===========================================================================
# RETURN STRUCTURE
# ===========================================================================

test_that("stepwise_reduce() returns a stepwise_reduction object", {
  expect_s3_class(LOG_SR, "stepwise_reduction")
})

test_that("stepwise_reduction object has all required fields", {
  expect_named(LOG_SR,
               c("final_model","original_model","step_table",
                 "dropped","kept","alpha","model_type"),
               ignore.order = TRUE)
})

test_that("model_type is 'Logistic' for glm", {
  expect_equal(LOG_SR$model_type, "Logistic")
})

test_that("model_type is 'Cox' for coxph", {
  expect_equal(COX_SR$model_type, "Cox")
})

test_that("alpha is stored correctly", {
  sr <- stepwise_reduce(LOG_FIT, alpha = 0.01, verbose = FALSE)
  expect_equal(sr$alpha, 0.01)
})

test_that("step_table has the expected columns", {
  expect_named(LOG_SR$step_table,
               c("step","term","df","AIC","dAIC",
                 "neg2LL","d_neg2LL","p","n_terms"),
               ignore.order = TRUE)
})

test_that("step_table rows equal number of dropped terms", {
  expect_equal(nrow(LOG_SR$step_table), length(LOG_SR$dropped))
})

test_that("step column is a contiguous integer sequence", {
  expect_equal(LOG_SR$step_table$step, seq_len(nrow(LOG_SR$step_table)))
})

test_that("n_terms decreases by one each step for numeric predictors", {
  n_start <- length(attr(terms(LOG_FIT), "term.labels"))
  expected <- seq(n_start - 1L, n_start - nrow(LOG_SR$step_table))
  expect_equal(LOG_SR$step_table$n_terms, expected)
})

test_that("dropped field matches step_table$term in order", {
  expect_equal(LOG_SR$dropped, LOG_SR$step_table$term)
})

test_that("final_model is a glm for logistic input", {
  expect_s3_class(LOG_SR$final_model, "glm")
})

test_that("final_model is a coxph for Cox input", {
  expect_s3_class(COX_SR$final_model, "coxph")
})

test_that("every dropped term had p >= alpha at the time it was removed", {
  # step_table$p records the LRT p-value that triggered each drop decision
  expect_true(all(LOG_SR$step_table$p >= 0.05))
})

test_that("every dropped Cox term had p >= alpha at the time it was removed", {
  expect_true(all(COX_SR$step_table$p >= 0.05))
})

# ===========================================================================
# LRT p-VALUES MATCH drop1() — CORE CORRECTNESS
# ===========================================================================

test_that("manual LRT p-values match drop1() for logistic (all numeric)", {
  d1_pv  <- drop1_pvals(LOG_FIT)
  man_pv <- manual_lrt_pvals(LOG_FIT, LOG_DF)
  common <- intersect(names(d1_pv), names(man_pv))
  expect_equal(man_pv[common], d1_pv[common], tolerance = 1e-8)
})

test_that("manual LRT p-values match drop1() for logistic (4-level factor)", {
  d1_pv  <- drop1_pvals(FAC_FIT)
  man_pv <- manual_lrt_pvals(FAC_FIT, FAC_DF)
  common <- intersect(names(d1_pv), names(man_pv))
  expect_equal(man_pv[common], d1_pv[common], tolerance = 1e-8)
})

test_that("manual LRT p-values match drop1() for Cox (all numeric)", {
  d1_pv  <- drop1_pvals(COX_FIT)
  man_pv <- manual_lrt_pvals(COX_FIT, COX_DF)
  common <- intersect(names(d1_pv), names(man_pv))
  expect_equal(man_pv[common], d1_pv[common], tolerance = 1e-8)
})

test_that("manual LRT p-values match drop1() for Cox (3-level factor)", {
  d1_pv  <- drop1_pvals(COXF_FIT)
  man_pv <- manual_lrt_pvals(COXF_FIT, COXF_DF)
  common <- intersect(names(d1_pv), names(man_pv))
  expect_equal(man_pv[common], d1_pv[common], tolerance = 1e-8)
})

# ===========================================================================
# FIRST-STEP DECISION MATCHES drop1()
# ===========================================================================

test_that("first term dropped matches drop1() for logistic (all numeric)", {
  expect_equal(LOG_SR$dropped[1L], drop1_first_drop(LOG_FIT))
})

test_that("first term dropped matches drop1() for logistic (4-level factor)", {
  expect_equal(FAC_SR$dropped[1L], drop1_first_drop(FAC_FIT))
})

test_that("first term dropped matches drop1() for Cox (all numeric)", {
  expect_equal(COX_SR$dropped[1L], drop1_first_drop(COX_FIT))
})

test_that("first term dropped matches drop1() for Cox (3-level factor)", {
  expect_equal(COXF_SR$dropped[1L], drop1_first_drop(COXF_FIT))
})

# ===========================================================================
# HIERARCHY: interactions
# ===========================================================================

test_that("main effects are not dropped before their interaction (step 1)", {
  # y ~ age * sex — only age:sex should be testable at step 1
  expect_equal(INT_SR$dropped[1L], "age:sex")
})

test_that("drop1() and stepwise_reduce agree on first testable term for interaction model", {
  expect_equal(INT_SR$dropped[1L], drop1_first_drop(INT_FIT))
})

test_that("after interaction is dropped, main effects become testable", {
  # INT_SR drops age:sex first; subsequent terms come from y ~ age + sex
  post_interaction <- INT_SR$dropped[-1L]
  expect_true(all(post_interaction %in% c("age", "sex")))
})

# ===========================================================================
# keep_vars
# ===========================================================================

test_that("keep_vars prevents a non-significant term from being dropped", {
  # bmi is not significant — protect it and confirm it survives
  sr_keep <- stepwise_reduce(LOG_FIT, alpha = 0.05,
                              keep_vars = "bmi", verbose = FALSE)
  expect_false("bmi" %in% sr_keep$dropped)
  expect_true("bmi" %in% attr(terms(sr_keep$final_model), "term.labels"))
})

test_that("keep_vars with Cox model protects a non-significant term", {
  sr_keep <- stepwise_reduce(COX_FIT, alpha = 0.05,
                              keep_vars = "noise", verbose = FALSE)
  expect_false("noise" %in% sr_keep$dropped)
})

test_that("multiple keep_vars are all protected", {
  sr_keep <- stepwise_reduce(LOG_FIT, alpha = 0.05,
                              keep_vars = c("bmi", "crp"), verbose = FALSE)
  expect_false("bmi"  %in% sr_keep$dropped)
  expect_false("crp"  %in% sr_keep$dropped)
})

# ===========================================================================
# NOTHING DROPPED WHEN ALL TERMS SIGNIFICANT
# ===========================================================================

test_that("no terms dropped when all are significant at alpha = 0.05", {
  # Fit model with only the significant predictor (age)
  fit_sig <- glm(y ~ age, data = LOG_DF, family = binomial)
  sr_sig  <- stepwise_reduce(fit_sig, alpha = 0.05, verbose = FALSE)
  expect_length(sr_sig$dropped, 0L)
  expect_equal(nrow(sr_sig$step_table), 0L)
})

test_that("no terms dropped when all Cox terms significant at alpha = 0.05", {
  # Lung data: sex and ph.ecog are clearly significant
  lung        <- survival::lung
  lung$status <- lung$status - 1L
  lung$sex    <- factor(lung$sex, labels = c("Male","Female"))
  lung_cc <- lung[complete.cases(lung[, c("time","status","sex","ph.ecog")]), ]
  fit_sig <- coxph(Surv(time, status) ~ sex + ph.ecog, data = lung_cc)
  sr_sig  <- stepwise_reduce(fit_sig, alpha = 0.05, verbose = FALSE)
  expect_length(sr_sig$dropped, 0L)
})

# ===========================================================================
# verbose = FALSE suppresses output
# ===========================================================================

test_that("verbose = FALSE produces no console output", {
  expect_silent(stepwise_reduce(LOG_FIT, alpha = 0.05, verbose = FALSE))
})

test_that("verbose = TRUE produces console output", {
  expect_output(stepwise_reduce(LOG_FIT, alpha = 0.05, verbose = TRUE))
})

# ===========================================================================
# S3 METHODS
# ===========================================================================

test_that("print.stepwise_reduction() runs without error", {
  expect_output(print(LOG_SR))
})

test_that("print.stepwise_reduction() returns the object invisibly", {
  out <- withVisible(print(LOG_SR))
  expect_false(out$visible)
  expect_identical(out$value, LOG_SR)
})

test_that("print output contains key header fields", {
  txt <- capture.output(print(LOG_SR))
  expect_true(any(grepl("Stepwise", txt)))
  expect_true(any(grepl("alpha",    txt)))
  expect_true(any(grepl("step",     txt, ignore.case = TRUE)))
})

test_that("summary.stepwise_reduction() runs without error", {
  expect_output(summary(LOG_SR))
})

test_that("summary.stepwise_reduction() returns the object invisibly", {
  out <- withVisible(summary(LOG_SR))
  expect_false(out$visible)
  expect_identical(out$value, LOG_SR)
})

test_that("summary output contains model comparison section", {
  txt <- capture.output(summary(LOG_SR))
  expect_true(any(grepl("Model Comparison|Original|Final", txt)))
  expect_true(any(grepl("AIC", txt)))
})

test_that("summary output contains LRT section", {
  txt <- capture.output(summary(LOG_SR))
  expect_true(any(grepl("LRT|chi2", txt, ignore.case = TRUE)))
})

test_that("print.stepwise_reduction() works for Cox result", {
  expect_output(print(COX_SR))
})

test_that("summary.stepwise_reduction() shows Cox-specific stats", {
  txt <- capture.output(summary(COX_SR))
  expect_true(any(grepl("Events|Cox|ties", txt, ignore.case = TRUE)))
})
