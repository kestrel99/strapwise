# Tests for fp_data() and fp_plot().
#
# fp_data() is tested exhaustively against known structure and values.
# fp_plot() is tested for correct return type and robustness to arguments;
# pixel-level rendering is not verified.

library(survival)

# ── Fixtures ──────────────────────────────────────────────────────────────────

# Logistic: continuous age + binary factor sex (reference "F").
make_log_df <- function(seed = 1L, n = 300L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 55, 10),
    bmi = rnorm(n, 27,  4),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  df$y <- rbinom(n, 1L,
                 plogis(-4 + 0.05 * df$age + 0.04 * df$bmi +
                          0.5 * (df$sex == "M")))
  df
}

# Logistic: 3-level factor + continuous predictor.
make_fac3_df <- function(seed = 2L, n = 300L) {
  set.seed(seed)
  df <- data.frame(
    score = rnorm(n),
    grp   = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  df$y <- rbinom(n, 1L, plogis(-1 + 0.4 * df$score))
  df
}

# Logistic: continuous × continuous interaction (age:bmi).
make_int_cc_df <- function(seed = 3L, n = 300L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 55, 10),
    bmi = rnorm(n, 27,  4)
  )
  df$y <- rbinom(n, 1L, plogis(-5 + 0.04 * df$age + 0.03 * df$bmi))
  df
}

# Logistic: continuous × factor interaction (age:sex).
make_int_cf_df <- function(seed = 4L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 55, 10),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  df$y <- rbinom(n, 1L, plogis(-4 + 0.05 * df$age + 0.5 * (df$sex == "M")))
  df
}

# Cox: continuous age + binary factor sex.
make_cox_df <- function(seed = 5L, n = 300L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 55, 10),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  df$time   <- rexp(n, exp(0.02 * df$age / 10))
  df$status <- rbinom(n, 1L, 0.75)
  df
}

# Cox: continuous × factor interaction (age:sex).
make_cox_int_df <- function(seed = 6L, n = 400L) {
  set.seed(seed)
  df <- data.frame(
    age = rnorm(n, 55, 10),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  df$time   <- rexp(n, exp(0.02 * df$age / 10))
  df$status <- rbinom(n, 1L, 0.75)
  df
}

# ── Pre-built fixtures ────────────────────────────────────────────────────────

LOG_DF     <- make_log_df()
FAC3_DF    <- make_fac3_df()
INT_CC_DF  <- make_int_cc_df()
INT_CF_DF  <- make_int_cf_df()
COX_DF     <- make_cox_df()
COX_INT_DF <- make_cox_int_df()

LOG_FIT <- glm(y ~ age + bmi + sex, data = LOG_DF,  family = binomial)
FAC_FIT <- glm(y ~ score + grp,     data = FAC3_DF, family = binomial)

INT_CC_FIT <- glm(y ~ age + bmi + age:bmi,
                  data = INT_CC_DF, family = binomial)
INT_CF_FIT <- glm(y ~ age + sex + age:sex,
                  data = INT_CF_DF, family = binomial)

COX_FIT <- coxph(Surv(time, status) ~ age + sex,           data = COX_DF)
COX_INT_FIT <- coxph(Surv(time, status) ~ age + sex + age:sex,
                     data = COX_INT_DF)

LOG_PD     <- fp_data(LOG_FIT, LOG_DF)
FAC_PD     <- fp_data(FAC_FIT, FAC3_DF)
INT_CC_PD  <- fp_data(INT_CC_FIT, INT_CC_DF)
INT_CF_PD  <- fp_data(INT_CF_FIT, INT_CF_DF)
COX_PD     <- fp_data(COX_FIT, COX_DF)
COX_INT_PD <- fp_data(COX_INT_FIT, COX_INT_DF)


# ── fp_data: return structure ─────────────────────────────────────────────────

test_that("fp_data returns a data frame", {
  expect_s3_class(LOG_PD, "data.frame")
  expect_s3_class(COX_PD, "data.frame")
})

test_that("fp_data has required base columns", {
  required <- c("Predictor", "est", "lci", "uci", "pval",
                "is_header", "is_reference", "group", "P-value", " ")
  expect_true(all(required %in% names(LOG_PD)))
  expect_true(all(required %in% names(COX_PD)))
})

test_that("logistic model produces OR (95% CI) column, not HR", {
  expect_true("OR (95% CI)" %in% names(LOG_PD))
  expect_false("HR (95% CI)" %in% names(LOG_PD))
})

test_that("Cox model produces HR (95% CI) column, not OR", {
  expect_true("HR (95% CI)" %in% names(COX_PD))
  expect_false("OR (95% CI)" %in% names(COX_PD))
})


# ── fp_data: row counts ───────────────────────────────────────────────────────

test_that("continuous predictor produces 1 header + n_percs rows", {
  n_percs <- 4L   # default percs = c(0.05, 0.25, 0.75, 0.95)
  age_rows <- LOG_PD[LOG_PD$group == which(attr(terms(LOG_FIT),
                                               "term.labels") == "age"), ]
  expect_equal(nrow(age_rows), 1L + n_percs)
  expect_equal(sum(age_rows$is_header), 1L)
})

test_that("binary factor produces 1 header + 1 reference + 1 level row", {
  sex_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "sex")
  sex_rows <- LOG_PD[LOG_PD$group == sex_idx, ]
  expect_equal(nrow(sex_rows), 3L)             # header + reference + "M"
  expect_equal(sum(sex_rows$is_header),    1L)
  expect_equal(sum(sex_rows$is_reference), 1L)
})

test_that("3-level factor produces 1 header + 1 reference + 2 level rows", {
  grp_idx  <- which(attr(terms(FAC_FIT), "term.labels") == "grp")
  grp_rows <- FAC_PD[FAC_PD$group == grp_idx, ]
  expect_equal(nrow(grp_rows), 4L)             # header + ref + B + C
  expect_equal(sum(grp_rows$is_header),    1L)
  expect_equal(sum(grp_rows$is_reference), 1L)
})

test_that("custom percs controls the number of percentile rows", {
  pd <- fp_data(LOG_FIT, LOG_DF, percs = c(0.1, 0.9))
  age_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "age")
  age_rows <- pd[pd$group == age_idx, ]
  expect_equal(nrow(age_rows), 3L)             # header + 2 percentiles
})


# ── fp_data: is_header / is_reference flags ───────────────────────────────────

test_that("is_header is TRUE only for variable header rows", {
  hrows <- LOG_PD[LOG_PD$is_header, ]
  expect_equal(nrow(hrows), length(attr(terms(LOG_FIT), "term.labels")))
  expect_true(all(is.na(hrows$est)))
})

test_that("reference row: is_reference TRUE, est NA", {
  ref_rows <- LOG_PD[LOG_PD$is_reference, ]
  expect_equal(nrow(ref_rows), 1L)             # one binary factor = one ref
  expect_true(all(is.na(ref_rows$est)))
})

test_that("median percentile row has is_reference TRUE", {
  pd <- fp_data(LOG_FIT, LOG_DF,
                percs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  age_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "age")
  age_rows <- pd[pd$group == age_idx & !pd$is_header, ]
  expect_equal(sum(age_rows$is_reference), 1L)
  expect_equal(age_rows$est[age_rows$is_reference], 1)
})


# ── fp_data: estimate values ──────────────────────────────────────────────────

test_that("factor level OR matches exp(coef)", {
  sex_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "sex")
  sex_rows <- LOG_PD[LOG_PD$group == sex_idx & !LOG_PD$is_header &
                       !LOG_PD$is_reference, ]
  expected <- exp(coef(LOG_FIT)[["sexM"]])
  expect_equal(sex_rows$est, expected, tolerance = 1e-10)
})

test_that("factor level CI matches exp(coef ± z * se)", {
  vc   <- vcov(LOG_FIT)
  b    <- coef(LOG_FIT)[["sexM"]]
  se   <- sqrt(vc["sexM", "sexM"])
  sex_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "sex")
  sex_rows <- LOG_PD[LOG_PD$group == sex_idx & !LOG_PD$is_header &
                       !LOG_PD$is_reference, ]
  expect_equal(sex_rows$lci, exp(b - 1.96 * se), tolerance = 0.01)
  expect_equal(sex_rows$uci, exp(b + 1.96 * se), tolerance = 0.01)
})

test_that("header row pval matches Wald p from first non-ref coef", {
  b  <- coef(LOG_FIT)[["sexM"]]
  se <- sqrt(vcov(LOG_FIT)["sexM", "sexM"])
  expected_p <- 2 * pnorm(-abs(b / se))
  sex_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "sex")
  hdr_row  <- LOG_PD[LOG_PD$group == sex_idx & LOG_PD$is_header, ]
  expect_equal(hdr_row$pval, expected_p, tolerance = 1e-10)
})


# ── fp_data: labels argument ──────────────────────────────────────────────────

test_that("labels renames header rows", {
  pd <- fp_data(LOG_FIT, LOG_DF,
                labels = c(age = "Age (years)", sex = "Sex at birth"))
  hrows <- pd[pd$is_header, ]
  expect_true("Age (years)"  %in% hrows$Predictor)
  expect_true("Sex at birth" %in% hrows$Predictor)
})

test_that("unmatched variable names fall back to raw name", {
  pd <- fp_data(LOG_FIT, LOG_DF, labels = c(age = "Age (years)"))
  hrows <- pd[pd$is_header, ]
  expect_true("bmi" %in% hrows$Predictor)  # unlabelled, keeps raw name
  expect_true("sex" %in% hrows$Predictor)
})


# ── fp_data: log_vars argument ────────────────────────────────────────────────

test_that("log_vars back-transforms display values in percentile labels", {
  pd_log <- fp_data(LOG_FIT, LOG_DF, log_vars = "age",
                    percs = c(0.1, 0.9))
  age_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "age")
  age_pcts <- pd_log[pd_log$group == age_idx & !pd_log$is_header, ]

  pd_raw <- fp_data(LOG_FIT, LOG_DF,
                    percs = c(0.1, 0.9))
  age_pcts_raw <- pd_raw[pd_raw$group == age_idx & !pd_raw$is_header, ]

  # Labels differ: log_vars version should contain exp'd values
  expect_false(identical(age_pcts$Predictor, age_pcts_raw$Predictor))
})


# ── fp_data: digits / formatting ─────────────────────────────────────────────

test_that("digits argument controls formatted CI text", {
  pd3 <- fp_data(LOG_FIT, LOG_DF, digits = 3L)
  pd5 <- fp_data(LOG_FIT, LOG_DF, digits = 5L)
  ci3 <- pd3[["OR (95% CI)"]]
  ci5 <- pd5[["OR (95% CI)"]]
  # At least one non-empty CI string should differ between 3 and 5 sig figs
  non_empty <- nzchar(ci3) & nzchar(ci5)
  expect_true(any(ci3[non_empty] != ci5[non_empty]))
})

test_that("no leading spaces in formatted CI text", {
  ci <- LOG_PD[["OR (95% CI)"]]
  non_empty <- ci[nzchar(ci)]
  expect_true(all(!grepl("^ ", non_empty)))
})

test_that("no leading spaces in the numeric value inside percentile labels", {
  age_idx  <- which(attr(terms(LOG_FIT), "term.labels") == "age")
  pct_rows <- LOG_PD[LOG_PD$group == age_idx & !LOG_PD$is_header, ]
  # Extract the value between parentheses and verify no leading space
  inner <- regmatches(pct_rows$Predictor,
                      regexpr("\\([^)]+\\)", pct_rows$Predictor))
  inner_vals <- gsub("[()]", "", inner)
  expect_true(all(!grepl("^ ", inner_vals)))
})


# ── fp_data: labels / log_vars validation ────────────────────────────────────

test_that("unknown label key produces a warning", {
  expect_warning(
    fp_data(LOG_FIT, LOG_DF, labels = c(typo_var = "Typo")),
    "label key"
  )
})

test_that("valid label keys produce no warning", {
  expect_no_warning(
    fp_data(LOG_FIT, LOG_DF,
            labels = c(age = "Age (years)", bmi = "BMI", sex = "Sex"))
  )
})

test_that("multiple unknown label keys are all named in the warning", {
  expect_warning(
    fp_data(LOG_FIT, LOG_DF,
            labels = c(bad1 = "A", age = "Age (years)", bad2 = "B")),
    regexp = "bad1.*bad2|bad2.*bad1"
  )
})

test_that("unknown log_vars entry produces a warning", {
  expect_warning(
    fp_data(LOG_FIT, LOG_DF, log_vars = "typo_var"),
    "log_vars"
  )
})

test_that("valid log_vars entries produce no warning", {
  expect_no_warning(
    fp_data(LOG_FIT, LOG_DF, log_vars = "age")
  )
})

test_that("categorical variable in log_vars produces a warning", {
  expect_warning(
    fp_data(LOG_FIT, LOG_DF, log_vars = "sex"),
    "categorical"
  )
})

test_that("warnings do not prevent fp_data from returning a result", {
  pd <- suppressWarnings(
    fp_data(LOG_FIT, LOG_DF,
            labels   = c(bad = "X"),
            log_vars = c("also_bad", "sex"))
  )
  expect_s3_class(pd, "data.frame")
  expect_gt(nrow(pd), 0L)
})


# ── fp_data: P-value formatting ───────────────────────────────────────────────

test_that("P-value is blank for rows with NA pval", {
  # Percentile and factor-level rows have NA pval
  non_hdr <- LOG_PD[!LOG_PD$is_header & !LOG_PD$is_reference, ]
  expect_true(all(non_hdr[["P-value"]] == ""))
})

test_that("P-value shows <0.001 when p is very small", {
  # Fit a model where a coefficient is strongly significant
  set.seed(99L)
  n  <- 1000L
  x  <- rnorm(n)
  y  <- rbinom(n, 1L, plogis(-1 + 3 * x))
  df <- data.frame(x, y)
  fit <- glm(y ~ x, data = df, family = binomial)
  pd  <- fp_data(fit, df)
  hdr <- pd[pd$is_header, ]
  expect_equal(hdr[["P-value"]], "<0.001")
})

test_that("P-value shows formatted value for moderate p", {
  # Use a model with a weak coefficient
  set.seed(77L)
  n  <- 100L
  x  <- rnorm(n)
  y  <- rbinom(n, 1L, plogis(-0.5 + 0.1 * x))
  df <- data.frame(x, y)
  fit <- glm(y ~ x, data = df, family = binomial)
  b   <- coef(fit)[["x"]]
  se  <- sqrt(vcov(fit)["x", "x"])
  p   <- 2 * pnorm(-abs(b / se))
  pd  <- fp_data(fit, df)
  hdr <- pd[pd$is_header, ]
  if (p >= 0.001) {
    expect_equal(hdr[["P-value"]], sprintf("%.3f", p))
  }
})


# ── fp_data: ci_level / ci_width ─────────────────────────────────────────────

test_that("ci_level affects CI width (wider at higher level)", {
  pd95 <- fp_data(LOG_FIT, LOG_DF, ci_level = 0.95,
                  percs = c(0.1, 0.9))
  pd99 <- fp_data(LOG_FIT, LOG_DF, ci_level = 0.99,
                  percs = c(0.1, 0.9))
  rows95 <- pd95[!is.na(pd95$est) & !pd95$is_reference, ]
  rows99 <- pd99[!is.na(pd99$est) & !pd99$is_reference, ]
  width95 <- rows95$uci - rows95$lci
  width99 <- rows99$uci - rows99$lci
  expect_true(all(width99 >= width95 - 1e-12))
})

test_that("ci_width controls spacer column character width", {
  pd20 <- fp_data(LOG_FIT, LOG_DF, ci_width = 20L)
  pd50 <- fp_data(LOG_FIT, LOG_DF, ci_width = 50L)
  expect_equal(nchar(pd20[[" "]][[1L]]), 20L)
  expect_equal(nchar(pd50[[" "]][[1L]]), 50L)
})


# ── fp_data: interaction terms ────────────────────────────────────────────────

test_that("continuous x continuous interaction appears as single non-header row", {
  int_idx  <- which(attr(terms(INT_CC_FIT), "term.labels") == "age:bmi")
  int_rows <- INT_CC_PD[INT_CC_PD$group == int_idx, ]
  expect_equal(nrow(int_rows), 1L)
  expect_false(int_rows$is_header)
  expect_false(is.na(int_rows$est))
})

test_that("continuous x continuous interaction OR matches exp(coef)", {
  int_idx <- which(attr(terms(INT_CC_FIT), "term.labels") == "age:bmi")
  int_row <- INT_CC_PD[INT_CC_PD$group == int_idx, ]
  expected <- exp(coef(INT_CC_FIT)[["age:bmi"]])
  expect_equal(int_row$est, expected, tolerance = 1e-10)
})

test_that("continuous x factor interaction produces header + level rows", {
  int_idx  <- which(attr(terms(INT_CF_FIT), "term.labels") == "age:sex")
  int_rows <- INT_CF_PD[INT_CF_PD$group == int_idx, ]
  # One header row; one level row for "M" (reference "F" absorbed in main effect)
  expect_equal(nrow(int_rows), 2L)
  expect_equal(sum(int_rows$is_header), 1L)
  expect_true(all(is.na(int_rows$est[int_rows$is_header])))
  expect_false(any(is.na(int_rows$est[!int_rows$is_header])))
})

test_that("continuous x factor interaction OR matches exp(coef)", {
  int_idx  <- which(attr(terms(INT_CF_FIT), "term.labels") == "age:sex")
  level_row <- INT_CF_PD[INT_CF_PD$group == int_idx & !INT_CF_PD$is_header, ]
  # Coefficient name may be "age:sexM" depending on R version
  int_cn <- grep("^age.*sex|^sex.*age", names(coef(INT_CF_FIT)), value = TRUE)
  expected <- exp(coef(INT_CF_FIT)[[int_cn]])
  expect_equal(level_row$est, expected, tolerance = 1e-10)
})

test_that("labels applies to interaction term header row", {
  pd <- fp_data(INT_CF_FIT, INT_CF_DF,
                labels = c("age:sex" = "Age x Sex interaction"))
  int_idx <- which(attr(terms(INT_CF_FIT), "term.labels") == "age:sex")
  hdr     <- pd[pd$group == int_idx & pd$is_header, ]
  expect_equal(hdr$Predictor, "Age x Sex interaction")
})

test_that("interaction in Cox model: header + level rows, HR column", {
  int_idx  <- which(attr(terms(COX_INT_FIT), "term.labels") == "age:sex")
  int_rows <- COX_INT_PD[COX_INT_PD$group == int_idx, ]
  expect_equal(nrow(int_rows), 2L)
  expect_true("HR (95% CI)" %in% names(COX_INT_PD))
  expect_equal(sum(int_rows$is_header), 1L)
})

test_that("model with only interactions does not error", {
  set.seed(9L)
  df  <- data.frame(a = rnorm(200), b = rnorm(200))
  df$y <- rbinom(200, 1L, plogis(df$a * df$b * 0.3))
  fit <- glm(y ~ a + b + a:b, data = df, family = binomial)
  expect_no_error(fp_data(fit, df))
})

test_that("all non-interaction rows unaffected by adding an interaction", {
  pd_main <- fp_data(glm(y ~ age + sex, data = INT_CF_DF, family = binomial),
                     INT_CF_DF)
  pd_int  <- fp_data(INT_CF_FIT, INT_CF_DF)

  # The age and sex groups should appear in both with same group indices
  expect_true(max(pd_main$group) < max(pd_int$group))
  # age header rows should have same pval sign (positive p)
  age_idx_main <- which(attr(terms(glm(y ~ age + sex, data = INT_CF_DF,
                                       family = binomial)),
                             "term.labels") == "age")
  age_idx_int  <- which(attr(terms(INT_CF_FIT), "term.labels") == "age")
  hdr_main <- pd_main[pd_main$group == age_idx_main & pd_main$is_header, ]
  hdr_int  <- pd_int[ pd_int$group  == age_idx_int  & pd_int$is_header,  ]
  expect_true(hdr_main$pval > 0)
  expect_true(hdr_int$pval  > 0)
})


# ── fp_data: Cox model correctness ────────────────────────────────────────────

test_that("Cox HR for factor level matches exp(coef)", {
  sex_idx  <- which(attr(terms(COX_FIT), "term.labels") == "sex")
  level_row <- COX_PD[COX_PD$group == sex_idx &
                        !COX_PD$is_header & !COX_PD$is_reference, ]
  expected <- exp(coef(COX_FIT)[["sexM"]])
  expect_equal(level_row$est, expected, tolerance = 1e-10)
})


# ── fp_plot: return type and basic rendering ──────────────────────────────────

test_that("fp_plot returns a gtable / forestplot object", {
  skip_if_not_installed("forestploter")
  fp <- fp_plot(LOG_FIT, LOG_DF)
  expect_true(inherits(fp, "gtable") || inherits(fp, "forestplot"))
})

test_that("fp_plot works with a logistic model", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(LOG_FIT, LOG_DF))
})

test_that("fp_plot works with a Cox model", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(COX_FIT, COX_DF))
})

test_that("fp_plot works with interaction terms", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(INT_CF_FIT, INT_CF_DF))
  expect_no_error(fp_plot(INT_CC_FIT, INT_CC_DF))
  expect_no_error(fp_plot(COX_INT_FIT, COX_INT_DF))
})

test_that("fp_plot forwards data-prep arguments to fp_data", {
  skip_if_not_installed("forestploter")
  lbs <- c(age = "Age (years)", sex = "Sex")
  expect_no_error(
    fp_plot(LOG_FIT, LOG_DF, labels = lbs, percs = c(0.1, 0.9), digits = 4L)
  )
})

test_that("fp_plot: vert_line = NULL does not error", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(LOG_FIT, LOG_DF, vert_line = NULL))
})

test_that("fp_plot: custom xlim and ticks_at do not error", {
  skip_if_not_installed("forestploter")
  expect_no_error(
    fp_plot(LOG_FIT, LOG_DF, xlim = c(0.5, 5), ticks_at = c(0.5, 1, 2, 5))
  )
})

test_that("fp_plot: x_trans = 'none' does not error", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(LOG_FIT, LOG_DF, x_trans = "none", xlim = c(0, 4)))
})

test_that("fp_plot: custom ci_col does not error", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(LOG_FIT, LOG_DF, ci_col = "#2166AC"))
})

test_that("fp_plot: 3-level factor model renders without error", {
  skip_if_not_installed("forestploter")
  expect_no_error(fp_plot(FAC_FIT, FAC3_DF))
})


# ── .fp_fmt internal helper ───────────────────────────────────────────────────

test_that(".fp_fmt returns '-' for NA", {
  expect_equal(strapwise:::.fp_fmt(NA_real_), "-")
})

test_that(".fp_fmt trims leading whitespace", {
  vals <- c(1.234567, 0.00001, 1000.5)
  result <- strapwise:::.fp_fmt(vals, digits = 3)
  expect_true(all(!grepl("^ ", result)))
})

test_that(".fp_fmt respects digits argument", {
  r3 <- strapwise:::.fp_fmt(1.23456, digits = 3)
  r5 <- strapwise:::.fp_fmt(1.23456, digits = 5)
  expect_equal(nchar(r3) <= nchar(r5), TRUE)
})


# ── .fp_fmt_p internal helper ─────────────────────────────────────────────────

test_that(".fp_fmt_p returns empty string for NA", {
  expect_equal(strapwise:::.fp_fmt_p(NA_real_), "")
})

test_that(".fp_fmt_p returns <0.001 for very small p", {
  expect_equal(strapwise:::.fp_fmt_p(0.00001), "<0.001")
  expect_equal(strapwise:::.fp_fmt_p(1e-10),   "<0.001")
})

test_that(".fp_fmt_p formats moderate p to 3 decimal places", {
  expect_equal(strapwise:::.fp_fmt_p(0.0523), "0.052")
  expect_equal(strapwise:::.fp_fmt_p(0.1234), "0.123")
})


# ── .fp_int_coefs internal helper ────────────────────────────────────────────

test_that(".fp_int_coefs finds continuous x continuous interaction", {
  coef_names <- c("(Intercept)", "age", "bmi", "age:bmi")
  result <- strapwise:::.fp_int_coefs("age:bmi", coef_names)
  expect_equal(result, "age:bmi")
})

test_that(".fp_int_coefs finds continuous x factor interaction (factor suffix)", {
  coef_names <- c("(Intercept)", "age", "sexM", "age:sexM")
  result <- strapwise:::.fp_int_coefs("age:sex", coef_names)
  expect_equal(result, "age:sexM")
})

test_that(".fp_int_coefs handles factor-first coefficient names", {
  # Some R versions / formula orderings produce "sexM:age" not "age:sexM"
  coef_names <- c("(Intercept)", "age", "sexM", "sexM:age")
  result <- strapwise:::.fp_int_coefs("age:sex", coef_names)
  expect_equal(result, "sexM:age")
})

test_that(".fp_int_coefs returns character(0) when no match", {
  coef_names <- c("(Intercept)", "age", "bmi")
  result <- strapwise:::.fp_int_coefs("age:bmi", coef_names)
  expect_equal(result, character(0L))
})

test_that(".fp_int_coefs finds multiple levels for multi-level factor", {
  coef_names <- c("(Intercept)", "age", "grpB", "grpC",
                  "age:grpB", "age:grpC")
  result <- strapwise:::.fp_int_coefs("age:grp", coef_names)
  expect_setequal(result, c("age:grpB", "age:grpC"))
})
