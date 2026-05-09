# =============================================================================
# strapwise — end-to-end logistic regression example
#
# Covers:
#   1. Standard logistic regression  (bootstrap_logistic / plot_logistic_curves
#      / compare_logistic_roc)
#   2. Logistic Emax (continuous predictor)  (fit_logistic_emax /
#      roc_logistic_emax / bootstrap_emax / plot_emax_curve)
# =============================================================================

library(strapwise)
library(ggplot2)


# =============================================================================
# PART 1 — Standard binomial GLM
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1  Simulate data
# -----------------------------------------------------------------------------
set.seed(101)
n   <- 400
age <- rnorm(n, mean = 55, sd = 12)
bmi <- rnorm(n, mean = 27, sd = 5)
sex <- factor(sample(c("F", "M"), n, replace = TRUE))

# True log-odds: baseline −5, age +0.05, bmi +0.07, male +0.6
lp  <- -5 + 0.05 * age + 0.07 * bmi + 0.6 * (sex == "M")
y   <- rbinom(n, 1, plogis(lp))

df  <- data.frame(y, age, bmi, sex)

# -----------------------------------------------------------------------------
# 1.2  Fit two nested models
# -----------------------------------------------------------------------------
fit_reduced <- glm(y ~ age + sex,       data = df, family = binomial)
fit_full    <- glm(y ~ age + bmi + sex, data = df, family = binomial)

summary(fit_full)

# -----------------------------------------------------------------------------
# 1.3  Bootstrap the full model
# -----------------------------------------------------------------------------
bl <- bootstrap_logistic(fit_full, n_boot = 500, seed = 42)
print(bl)
summary(bl)

# -----------------------------------------------------------------------------
# 1.4  Response-curve plots (one per predictor)
# -----------------------------------------------------------------------------
plots <- plot_logistic_curves(bl)

plots[["age"]]   # continuous predictor
plots[["bmi"]]
plots[["sex"]]   # categorical predictor (dot + CI bars)

# Colour raw jitter points by sex; split observed-proportion overlay by group
plots_grp <- plot_logistic_curves(
  bl,
  colour_by      = "sex",
  obs_by_group   = TRUE,
  obs_groups     = 4L,
  colour_by_label = "Sex"
)
plots_grp[["age"]]
plots_grp[["bmi"]]

# Combine all panels into one figure (requires patchwork)
if (requireNamespace("patchwork", quietly = TRUE)) {
  plot_logistic_curves(bl, combine = TRUE)
}

# -----------------------------------------------------------------------------
# 1.5  ROC comparison
# -----------------------------------------------------------------------------
compare_logistic_roc(
  list("Age + sex"    = fit_reduced,
       "Full model"   = fit_full),
  legend_title = "Model"
)


# =============================================================================
# PART 2 — Logistic Emax model (continuous predictor)
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1  Simulate Emax data (dose as predictor) with a linear covariate
# -----------------------------------------------------------------------------
set.seed(202)
m    <- 300
dose <- sample(c(0, 5, 10, 25, 50, 100), m, replace = TRUE)
age2 <- rnorm(m, mean = 45, sd = 10)

# True model: E0 = −2, Emax = 3, EC50 = 15, beta_age = 0.03
eta2 <- -2 + 3 * dose / (15 + dose) + 0.03 * (age2 - mean(age2))
y2   <- rbinom(m, 1, plogis(eta2))

# -----------------------------------------------------------------------------
# 2.2  Fit models: predictor only, predictor + linear covariate
# -----------------------------------------------------------------------------
fit_emax_base <- fit_logistic_emax(y2, dose)
print(fit_emax_base)

fit_emax_cov  <- fit_logistic_emax(
  y2, dose,
  linear_covs = cbind(age_c = age2 - mean(age2))
)
print(fit_emax_cov)

# Delta-method confidence intervals for a new subject (age = 50)
nd_age  <- cbind(age_c = 50 - mean(age2))
predict(fit_emax_cov, newlinear = nd_age, type = "ci")

# -----------------------------------------------------------------------------
# 2.3  ROC / discrimination
# -----------------------------------------------------------------------------
roc_base <- roc_logistic_emax(fit_emax_base)
print(roc_base)

roc_cov  <- roc_logistic_emax(fit_emax_cov)
print(roc_cov)

# -----------------------------------------------------------------------------
# 2.4  Bootstrap
# -----------------------------------------------------------------------------
be_base <- bootstrap_emax(fit_emax_base, n_boot = 500, seed = 42)
print(be_base)
summary(be_base)

be_cov  <- bootstrap_emax(fit_emax_cov, n_boot = 500, seed = 42)
print(be_cov)
summary(be_cov)

# -----------------------------------------------------------------------------
# 2.5  Response curve plot (base model)
# -----------------------------------------------------------------------------
plot_emax_curve(be_base)

# With observed proportions and jitter; supply a domain-specific x label
plot_emax_curve(
  be_base,
  x_range    = c(0, 120),
  obs_groups = 5L,
  raw_data   = TRUE,
  x_label    = "Dose (mg)",
  y_label    = "P(Response)"
)

# Covariate model: curve evaluated at age = 45 (mean) and age = 65
nd_mean <- cbind(age_c = 0)
nd_old  <- cbind(age_c = 20)

p_mean <- plot_emax_curve(
  be_cov, newlinear = matrix(rep(0,  200), ncol = 1),
  x_label = "Dose (mg)", y_label = "P(Response)", raw_data = FALSE
)
p_old  <- plot_emax_curve(
  be_cov, newlinear = matrix(rep(20, 200), ncol = 1),
  x_label = "Dose (mg)", y_label = "P(Response)", raw_data = FALSE
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  (p_mean + ggplot2::ggtitle("Age = 45")) +
  (p_old  + ggplot2::ggtitle("Age = 65"))
}
