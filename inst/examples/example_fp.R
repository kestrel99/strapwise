# =============================================================================
# Forest plot helpers: end-to-end example
#
# Covers:
#   fp_data()   — build a tidy data frame from a glm or coxph model
#                 (useful for inspecting or modifying data before plotting)
#   fp_plot()   — build and render a forest plot in one step
#
# Two worked examples:
#   PART 1 — Logistic regression (binomial GLM)
#   PART 2 — Cox proportional hazards
# =============================================================================

library(strapwise)
library(survival)


# =============================================================================
# PART 1 — Logistic regression
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1  Simulate data
# -----------------------------------------------------------------------------
# Mix of predictor types:
#   age_log  — age stored on the log scale (continuous, log_vars demo)
#   bmi      — BMI, continuous, natural scale
#   ecog     — ECOG PS 0/1, binary factor
# True log-OR: intercept −6, log(age) +1.8, bmi +0.06, ECOG≥1 +0.8

set.seed(42)
n       <- 500
age     <- rnorm(n, mean = 60, sd = 12)
age_log <- log(pmax(age, 1))
bmi     <- rnorm(n, mean = 27, sd = 5)
ecog    <- factor(sample(c("0", "1+"), n, replace = TRUE, prob = c(0.6, 0.4)))

lp <- -6 + 1.8 * age_log + 0.06 * bmi + 0.8 * (ecog == "1+")
y  <- rbinom(n, 1L, plogis(lp))

logistic_dat <- data.frame(y, age_log, bmi, ecog)

# -----------------------------------------------------------------------------
# 1.2  Fit model
# -----------------------------------------------------------------------------
logistic_fit <- glm(
  y ~ age_log + bmi + ecog,
  data   = logistic_dat,
  family = binomial
)
summary(logistic_fit)

# -----------------------------------------------------------------------------
# 1.3  (Optional) inspect the intermediate data frame
# -----------------------------------------------------------------------------
# fp_data() is available separately when you need to examine or modify the
# plot data before rendering.  `log_vars` back-transforms stored log-scale
# values so percentile labels show original units (years, not log-years).

logistic_labels <- c(
  age_log = "Age (years)",
  bmi     = "BMI (kg/m²)",
  ecog    = "Baseline ECOG PS"
)

logistic_pd <- fp_data(
  model    = logistic_fit,
  data     = logistic_dat,
  labels   = logistic_labels,
  log_vars = "age_log",
  percs    = c(0.1, 0.25, 0.75, 0.9)
)

print(logistic_pd[, c("Predictor", "OR (95% CI)", "P-value", "is_header")])

# -----------------------------------------------------------------------------
# 1.4  Render the forest plot
# -----------------------------------------------------------------------------
# fp_plot() calls fp_data() internally — pass all arguments in one step.

fp_plot(
  logistic_fit,
  logistic_dat,
  labels   = logistic_labels,
  log_vars = "age_log",
  percs    = c(0.1, 0.25, 0.75, 0.9)
)

# Tighter x limits and custom axis ticks
fp_plot(
  logistic_fit, logistic_dat,
  labels   = logistic_labels,
  log_vars = "age_log",
  xlim     = c(0.5, 5),
  ticks_at = c(0.5, 1, 2, 5)
)

# Wider margin lines and a different CI colour
fp_plot(
  logistic_fit, logistic_dat,
  labels    = logistic_labels,
  log_vars  = "age_log",
  vert_line = c(0.5, 2),
  ci_col    = "#2166AC"
)

# Linear (non-log) x-axis — useful when all ORs are close to 1
fp_plot(
  logistic_fit, logistic_dat,
  labels   = logistic_labels,
  log_vars = "age_log",
  x_trans  = "none",
  xlim     = c(0, 4)
)


# =============================================================================
# PART 2 — Cox proportional hazards
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1  Data: NCCTG lung cancer dataset
# -----------------------------------------------------------------------------
lung <- survival::lung

lung$status   <- lung$status - 1L
lung$sex      <- factor(lung$sex, levels = 1:2, labels = c("Male", "Female"))
lung$ecog_cat <- factor(
  ifelse(lung$ph.ecog >= 2L, "2+", as.character(lung$ph.ecog)),
  levels = c("0", "1", "2+")
)
lung$meal_log <- log(pmax(lung$meal.cal, 1))

cox_dat <- lung[complete.cases(lung[, c("time", "status",
                                        "age", "sex",
                                        "ecog_cat", "meal_log")]), ]
nrow(cox_dat)

# -----------------------------------------------------------------------------
# 2.2  Fit Cox model
# -----------------------------------------------------------------------------
cox_fit <- coxph(
  Surv(time, status) ~ age + sex + ecog_cat + meal_log,
  data = cox_dat
)
summary(cox_fit)

# -----------------------------------------------------------------------------
# 2.3  Render the forest plot in one step
# -----------------------------------------------------------------------------
cox_labels <- c(
  age      = "Age (years)",
  sex      = "Sex",
  ecog_cat = "Baseline ECOG PS",
  meal_log = "Meal calories (kcal)"
)

fp_plot(
  cox_fit, cox_dat,
  labels   = cox_labels,
  log_vars = "meal_log",
  percs    = c(0.1, 0.25, 0.75, 0.9)
)

# Narrower x range; extra tick marks at standard HR reference points
fp_plot(
  cox_fit, cox_dat,
  labels   = cox_labels,
  log_vars = "meal_log",
  xlim     = c(0.3, 3),
  ticks_at = c(0.3, 0.5, 1, 2, 3)
)

# Suppress the ±20 % margin lines
fp_plot(
  cox_fit, cox_dat,
  labels    = cox_labels,
  log_vars  = "meal_log",
  vert_line = NULL,
  ci_col    = "#1A9641"
)


# =============================================================================
# 3.  Saving the output
# =============================================================================

fp <- fp_plot(cox_fit, cox_dat, labels = cox_labels, log_vars = "meal_log")

png("cox_forest.png", width = 9, height = 5, units = "in", res = 300)
grid::grid.draw(fp)
dev.off()

pdf("cox_forest.pdf", width = 9, height = 5)
grid::grid.draw(fp)
dev.off()
