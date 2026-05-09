# =============================================================================
# strapwise — end-to-end Cox proportional hazards example
#
# Covers:
#   bootstrap_cox()      — case-resampling bootstrap for a coxph fit
#   print / summary      — inspect results
#   plot_cox_forest()    — forest plot of hazard ratios with bootstrap CIs
#   plot_cox_survival()  — adjusted survival curves with bootstrap CI ribbons
# =============================================================================

library(strapwise)
library(survival)
library(ggplot2)


# =============================================================================
# 1.  Data
# =============================================================================
# Use the NCCTG lung cancer dataset bundled with the survival package.
# status is recoded 0/1 (survival convention: 0 = censored, 1 = event).

data(lung, package = "survival")
lung$status <- lung$status - 1L          # recode: 1→0 (censored), 2→1 (dead)
lung$sex    <- factor(lung$sex,           # 1 = male, 2 = female
                      levels = 1:2,
                      labels = c("Male", "Female"))

# Drop rows missing any predictor we plan to use
lung_cc <- lung[complete.cases(lung[, c("time", "status",
                                        "age", "sex", "ph.ecog")]), ]
nrow(lung_cc)   # 227 complete cases


# =============================================================================
# 2.  Fit Cox model
# =============================================================================
fit <- coxph(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung_cc
)
summary(fit)


# =============================================================================
# 3.  Bootstrap
# =============================================================================
bc <- bootstrap_cox(fit, n_boot = 1000, seed = 42)

print(bc)     # compact header + original log-HRs
summary(bc)   # HR, SE, bias, 95 % percentile CI on both scales


# =============================================================================
# 4.  Forest plot
# =============================================================================

# Default: log scale x-axis, HR = 1 reference line
plot_cox_forest(bc)

# Linear scale (HRs plotted on arithmetic axis)
plot_cox_forest(bc, log_scale = FALSE)

# Wider points and lines, custom colour
plot_cox_forest(
  bc,
  point_size = 4,
  line_size  = 1.0,
  palette    = c("#117733")
)


# =============================================================================
# 5.  Adjusted survival curves
# =============================================================================

# Define covariate profiles to compare
newdata <- data.frame(
  age     = c(50,  70,  50,  70),
  sex     = factor(c("Male", "Male", "Female", "Female"),
                   levels = c("Male", "Female")),
  ph.ecog = c(0,   2,   0,   2)
)

labels <- c(
  "Age 50, Male, ECOG 0",
  "Age 70, Male, ECOG 2",
  "Age 50, Female, ECOG 0",
  "Age 70, Female, ECOG 2"
)

# ---- 5a. Survival probability S(t)
plot_cox_survival(
  bc,
  newdata      = newdata,
  curve_labels = labels,
  x_label      = "Days",
  y_label      = "Survival Probability"
)

# ---- 5b. Event probability 1 − S(t)
plot_cox_survival(
  bc,
  newdata      = newdata,
  curve_labels = labels,
  fun          = "event",
  x_label      = "Days"
)

# ---- 5c. Cumulative hazard H(t) = −log S(t)
plot_cox_survival(
  bc,
  newdata      = newdata,
  curve_labels = labels,
  fun          = "cumhaz",
  x_label      = "Days"
)

# ---- 5d. Evaluate on a regular time grid (0 to 1000 days, every 10 days)
t_grid <- seq(0, 1000, by = 10)

plot_cox_survival(
  bc,
  newdata      = newdata[1:2, ],
  curve_labels = labels[1:2],
  times        = t_grid,
  x_label      = "Days"
)

# ---- 5e. 90 % CI, thinner ribbons
plot_cox_survival(
  bc,
  newdata      = newdata,
  curve_labels = labels,
  ribbon_alpha = 0.10,
  x_label      = "Days"
)


# =============================================================================
# 6.  Combine forest + survival with patchwork
# =============================================================================
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)

  p_forest <- plot_cox_forest(bc)
  p_surv   <- plot_cox_survival(
    bc,
    newdata      = newdata,
    curve_labels = labels,
    x_label      = "Days"
  )

  p_forest + p_surv + plot_layout(widths = c(1, 2))
}
