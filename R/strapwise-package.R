#' strapwise: Bootstrap Inference and Response Curves for Logistic Regression
#'
#' Provides case-resampling bootstrap for fitted binomial GLMs and a companion
#' plotting function that generates marginal response curves for all predictors
#' in the model.
#'
#' The main entry points are:
#'
#' - [bootstrap_logistic()]: bootstrap a fitted `glm` object.
#' - [compare_logistic_roc()]: compare multiple models with overlaid ROC curves.
#' - [plot_logistic_curves()]: generate response-curve plots from the result.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats AIC BIC coef formula glm model.frame plogis pnorm qnorm
#'   quantile rbinom rnorm setNames
## usethis namespace: end
NULL

## Suppress R CMD CHECK NOTEs for ggplot2 aesthetic names used inside aes()
utils::globalVariables(
  c(
    "x", "y", "prob", "ci_lo", "ci_hi", "obs", "grp", "D", "M", "model",
    "fpr", "sensitivity", "threshold", "specificity",
    "term", "hr", "time", "surv", "curve"
  )
)
