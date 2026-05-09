# strapwise

Bootstrap inference, response-curve plotting, and forest plots for logistic
regression and Cox proportional hazards models.

## Overview

`strapwise` provides a toolkit for fitting, bootstrapping, and visualising
binary-outcome and time-to-event regression models.

### Logistic regression (`glm(..., family = binomial)`)

| Function | Description |
|---|---|
| `bootstrap_logistic()` | Case-resampling bootstrap; returns coefficient draws and percentile CIs |
| `plot_logistic_curves()` | Marginal response-curve plots for each predictor with bootstrap CI ribbon |
| `compare_logistic_roc()` | Overlaid ROC curves from multiple fitted GLMs |
| `stepwise_reduce()` | Backwards stepwise reduction by likelihood ratio test |

### Logistic Emax models

| Function | Description |
|---|---|
| `fit_logistic_emax()` | Fit a logistic Emax dose-response model via maximum likelihood |
| `bootstrap_emax()` | Case-resampling bootstrap for a `logistic_emax` object |
| `plot_emax_curve()` | Response-probability curve with bootstrap CI ribbon |
| `roc_logistic_emax()` | ROC curve with AUC and Youden-optimal threshold |

### Cox proportional hazards models (`survival::coxph()`)

| Function | Description |
|---|---|
| `bootstrap_cox()` | Case-resampling bootstrap; returns log-HR draws and percentile CIs |
| `plot_cox_forest()` | Forest plot of HR point estimates and bootstrap CI bars |
| `plot_cox_survival()` | Adjusted survival curves with bootstrap CI ribbons |
| `stepwise_reduce()` | Also supports `coxph` objects |

### Forest plots (`glm` and `coxph`)

| Function | Description |
|---|---|
| `fp_plot()` | One-step forest plot from a fitted `glm` or `coxph` object |
| `fp_data()` | Extract the intermediate data frame before plotting |

## Install

```r
# from GitHub
# install.packages("pak")
pak::pak("kestrel99/strapwise")

# or with remotes
# install.packages("remotes")
remotes::install_github("kestrel99/strapwise")
```

## Examples

### Logistic regression bootstrap

```r
library(strapwise)

df <- data.frame(
  y   = rbinom(200, 1, 0.4),
  age = rnorm(200, 50, 10),
  sex = factor(sample(c("F", "M"), 200, replace = TRUE))
)

fit <- glm(y ~ age + sex, data = df, family = binomial)
bl  <- bootstrap_logistic(fit, n_boot = 500, seed = 1)

summary(bl)
plot_logistic_curves(bl)[["age"]]
```

### Forest plot (logistic)

```r
fp_plot(fit, df,
  labels = c(age = "Age (years)", sex = "Sex"))
```

### Cox model bootstrap and survival curves

```r
library(survival)

lung <- survival::lung
lung$status <- lung$status - 1L
fit_cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)

bc <- bootstrap_cox(fit_cox, n_boot = 500, seed = 1)
summary(bc)

plot_cox_forest(bc)

nd <- data.frame(age = c(50, 70), sex = c(1, 1), ph.ecog = c(0, 2))
plot_cox_survival(bc, newdata = nd,
                  curve_labels = c("Age 50, ECOG 0", "Age 70, ECOG 2"))
```

### Forest plot (Cox)

```r
fp_plot(fit_cox, lung,
  labels = c(age = "Age (years)", sex = "Sex", ph.ecog = "ECOG PS"))
```

### Logistic Emax dose-response

```r
set.seed(42)
dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = 200))
y    <- rbinom(200, 1, plogis(qlogis(0.10) +
         (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)))

fit_emax <- fit_logistic_emax(y, dose)
be       <- bootstrap_emax(fit_emax, n_boot = 500, seed = 1)

plot_emax_curve(be)
roc_logistic_emax(fit_emax)
```
