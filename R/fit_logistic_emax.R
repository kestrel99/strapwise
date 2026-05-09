#' Fit a logistic Emax model
#'
#' Fits a logistic Emax model to binary outcome data via maximum likelihood,
#' with optional additive linear covariates and covariate effects on the Emax
#' and EC50 parameters.
#'
#' The model for observation \eqn{i} is
#' \deqn{
#'   \eta_i = E_0
#'     + \bigl(E_{max} + \mathbf{z}_{emax,i}^\top\boldsymbol{\gamma}\bigr)
#'       \frac{x_i}{e^{\ell_i} + x_i}
#'     + \mathbf{z}_{lin,i}^\top\boldsymbol{\beta},
#' }
#' where \eqn{\ell_i = \log(EC_{50}) + \mathbf{z}_{ec50,i}^\top\boldsymbol{\delta}}
#' keeps \eqn{EC_{50,i} > 0}, and
#' \eqn{P(Y_i=1) = \text{logistic}(\eta_i)}.
#' When all covariate arguments are `NULL` the model reduces to the standard
#' three-parameter logistic Emax formula.
#'
#' @section Covariate matrices:
#' `linear_covs`, `emax_covs`, and `ec50_covs` must be numeric matrices (or
#' objects coercible via [as.matrix()]) with `length(y)` rows.  Expand factors
#' in advance with [model.matrix()].  Column names label the coefficients;
#' unnamed columns receive labels `"V1"`, `"V2"`, etc.
#'
#' @param y Numeric binary outcome vector (0 or 1).
#' @param x Numeric predictor vector (e.g. dose, concentration, time), same
#'   length as `y`.
#' @param linear_covs Optional numeric matrix of covariates that enter the
#'   linear predictor additively (\eqn{\mathbf{z}_{lin}}). `NULL` omits.
#' @param emax_covs Optional numeric matrix of covariates that modify the Emax
#'   parameter additively (\eqn{\mathbf{z}_{emax}}). `NULL` omits.
#' @param ec50_covs Optional numeric matrix of covariates that modify
#'   \eqn{\log(EC_{50})} additively (\eqn{\mathbf{z}_{ec50}}); a coefficient
#'   of \eqn{\delta} multiplies EC50 by \eqn{e^\delta} per unit increase.
#'   `NULL` omits.
#' @param start Starting values. One of:
#'   \itemize{
#'     \item `NULL` (default): auto-derived from a simple logistic regression;
#'       covariate effects start at 0.
#'     \item A length-3 vector `c(e0, emax, ec50)` on the natural scale:
#'       base parameters; covariate effects start at 0.
#'     \item A full-length vector of length
#'       `3 + ncol(linear_covs) + ncol(emax_covs) + ncol(ec50_covs)` ordered
#'       `(e0, emax, ec50, beta..., gamma..., delta...)`, all on the natural
#'       scale (ec50 is log-transformed internally).
#'   }
#' @param control A list of control parameters passed to [optim()]. Default
#'   `list(maxit = 1000)`.
#' @param hessian Logical; if `TRUE` (default) the Hessian is computed and used
#'   to derive the variance-covariance matrix and standard errors.
#' @param ... Additional arguments passed to [optim()].
#'
#' @return An S3 object of class `"logistic_emax"`, a named list with:
#'   \describe{
#'     \item{`coefficients`}{Named numeric vector on the reporting scale:
#'       `c(e0, emax, ec50, linear.*, emax.*, log_ec50.*)`.  `ec50` is
#'       back-transformed; EC50 covariate effects (`log_ec50.*`) remain on the
#'       log scale.}
#'     \item{`std_err`}{Standard errors matching `coefficients` (delta method
#'       applied to `ec50`).}
#'     \item{`vcov_theta`}{Variance-covariance matrix in the optimisation
#'       parameterisation `(e0, emax0, log_ec50_0, beta..., gamma...,
#'       delta...)`, or `NULL` when unavailable.}
#'     \item{`theta_hat`}{Raw optimised parameter vector (same parameterisation
#'       as `vcov_theta`).}
#'     \item{`fitted.values`}{Fitted probabilities for the observed data.}
#'     \item{`linear.predictors`}{Fitted values on the logit scale.}
#'     \item{`residuals`}{Raw residuals `y - fitted.values`.}
#'     \item{`loglik`}{Log-likelihood at the MLE.}
#'     \item{`AIC`}{Akaike information criterion.}
#'     \item{`converged`}{`TRUE` when [optim()] converged.}
#'     \item{`n_par`}{Total number of parameters.}
#'     \item{`x`, `y`}{Complete-case predictor and outcome vectors.}
#'     \item{`Z_lin`, `Z_emax`, `Z_ec50`}{Covariate matrices used in fitting
#'       (0-column matrices when the corresponding argument was `NULL`).}
#'     \item{`idx`}{Named list of integer indices into `theta_hat` for each
#'       parameter group, used internally by [predict.logistic_emax()].}
#'     \item{`call`}{The matched function call.}
#'   }
#'
#' @seealso [predict.logistic_emax()], [print.logistic_emax()],
#'   [roc_logistic_emax()]
#'
#' @examples
#' set.seed(42)
#' n    <- 200
#' dose <- sort(rep(c(0, 5, 10, 25, 50, 100), length.out = n))
#' eta  <- qlogis(0.10) + (qlogis(0.75) - qlogis(0.10)) * dose / (10 + dose)
#' y    <- rbinom(n, 1, plogis(eta))
#'
#' # Basic model
#' fit <- fit_logistic_emax(y, dose)
#' print(fit)
#'
#' # Linear covariate on the link scale
#' wt  <- rnorm(n, 70, 10)
#' fit2 <- fit_logistic_emax(y, dose, linear_covs = cbind(weight = wt))
#' print(fit2)
#'
#' # Covariate modifying Emax
#' grp <- rbinom(n, 1, 0.5)
#' fit3 <- fit_logistic_emax(y, dose, emax_covs = cbind(group = grp))
#' print(fit3)
#'
#' @export
fit_logistic_emax <- function(
  y,
  x,
  linear_covs = NULL,
  emax_covs   = NULL,
  ec50_covs   = NULL,
  start       = NULL,
  control     = list(maxit = 1000),
  hessian     = TRUE,
  ...
) {
  if (length(y) != length(x)) stop("y and x must have the same length.")
  if (!all(y %in% c(0, 1)))   stop("y must be binary (0/1).")

  # ---- Covariate matrices (original length, before subsetting) ---------------
  to_mat <- function(z, n_orig, label) {
    if (is.null(z)) return(matrix(0, nrow = n_orig, ncol = 0))
    m <- as.matrix(z)
    if (!is.numeric(m))    stop(label, " must be numeric.")
    if (nrow(m) != n_orig) stop(label, " must have ", n_orig, " rows.")
    if (ncol(m) > 0 && is.null(colnames(m)))
      colnames(m) <- paste0("V", seq_len(ncol(m)))
    m
  }
  n_orig     <- length(y)
  Z_lin_raw  <- to_mat(linear_covs, n_orig, "linear_covs")
  Z_emax_raw <- to_mat(emax_covs,   n_orig, "emax_covs")
  Z_ec50_raw <- to_mat(ec50_covs,   n_orig, "ec50_covs")

  # ---- Complete-case subsetting ----------------------------------------------
  ok <- stats::complete.cases(y, x)
  if (ncol(Z_lin_raw)  > 0) ok <- ok & stats::complete.cases(Z_lin_raw)
  if (ncol(Z_emax_raw) > 0) ok <- ok & stats::complete.cases(Z_emax_raw)
  if (ncol(Z_ec50_raw) > 0) ok <- ok & stats::complete.cases(Z_ec50_raw)

  y      <- y[ok]; x <- x[ok]; n <- length(y)
  Z_lin  <- Z_lin_raw[ok,  , drop = FALSE]
  Z_emax <- Z_emax_raw[ok, , drop = FALSE]
  Z_ec50 <- Z_ec50_raw[ok, , drop = FALSE]

  p_lin  <- ncol(Z_lin)
  p_emax <- ncol(Z_emax)
  p_ec50 <- ncol(Z_ec50)
  n_par  <- 3L + p_lin + p_emax + p_ec50

  # ---- Parameter indices (1-based into theta) --------------------------------
  idx_e0       <- 1L
  idx_emax0    <- 2L
  idx_log_ec50 <- 3L
  idx_lin      <- if (p_lin  > 0) seq.int(4L,             3L + p_lin)                    else integer(0)
  idx_emax_cov <- if (p_emax > 0) seq.int(4L + p_lin,     3L + p_lin + p_emax)           else integer(0)
  idx_ec50_cov <- if (p_ec50 > 0) seq.int(4L + p_lin + p_emax, n_par)                    else integer(0)

  idx <- list(
    e0       = idx_e0,
    emax0    = idx_emax0,
    log_ec50 = idx_log_ec50,
    lin      = idx_lin,
    emax_cov = idx_emax_cov,
    ec50_cov = idx_ec50_cov
  )

  # ---- Parameter names (optimisation scale) ----------------------------------
  theta_names <- c(
    "e0", "emax", "log_ec50",
    if (p_lin  > 0) paste0("linear.",   colnames(Z_lin))  else character(0),
    if (p_emax > 0) paste0("emax.",     colnames(Z_emax)) else character(0),
    if (p_ec50 > 0) paste0("log_ec50.", colnames(Z_ec50)) else character(0)
  )

  # ---- Negative log-likelihood -----------------------------------------------
  nll <- function(theta) {
    emax_i <- theta[idx_emax0] +
      if (p_emax > 0) drop(Z_emax %*% theta[idx_emax_cov]) else 0
    ec50_i <- exp(theta[idx_log_ec50] +
      if (p_ec50 > 0) drop(Z_ec50 %*% theta[idx_ec50_cov]) else 0)
    eta <- theta[idx_e0] + emax_i * x / (ec50_i + x) +
      if (p_lin > 0) drop(Z_lin %*% theta[idx_lin]) else 0
    -sum(
      y       * stats::plogis(eta, log.p = TRUE) +
      (1 - y) * stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
    )
  }

  # ---- Starting values -------------------------------------------------------
  if (is.null(start)) {
    fit_lin <- try(
      suppressWarnings(stats::glm(y ~ x, family = stats::binomial)),
      silent = TRUE
    )
    if (!inherits(fit_lin, "try-error") && all(is.finite(stats::coef(fit_lin)))) {
      b0         <- stats::coef(fit_lin)[[1]]
      b1         <- stats::coef(fit_lin)[[2]]
      emax_start <- b1 * diff(range(x))
      if (abs(emax_start) < 1e-4) emax_start <- 1
      e0_start   <- b0
    } else {
      e0_start   <- stats::qlogis(mean(pmin(pmax(y, 0.01), 0.99)))
      emax_start <- 1
    }
    x_pos      <- unique(x[x > 0 & is.finite(x)])
    ec50_start <- if (length(x_pos) > 0) stats::median(x_pos) else 1
    start_full <- c(e0_start, emax_start, log(ec50_start), rep(0, n_par - 3L))

  } else if (length(start) == 3L) {
    start_full <- c(start[[1]], start[[2]], log(start[[3]]), rep(0, n_par - 3L))

  } else if (length(start) == n_par) {
    start_full              <- as.numeric(start)
    start_full[idx_log_ec50] <- log(start_full[idx_log_ec50])

  } else {
    stop(sprintf(
      "`start` must have length 3 or %d (= total number of parameters).", n_par
    ))
  }
  names(start_full) <- theta_names

  # ---- Optimise --------------------------------------------------------------
  opt <- stats::optim(
    par     = start_full,
    fn      = nll,
    method  = "BFGS",
    control = control,
    hessian = hessian,
    ...
  )
  if (opt$convergence != 0)
    warning("Optimisation did not converge. Code: ", opt$convergence)

  th         <- setNames(opt$par, theta_names)
  ec50_est   <- exp(th[idx_log_ec50])

  # ---- Variance-covariance ---------------------------------------------------
  vcov_theta <- NULL
  if (hessian && all(is.finite(opt$hessian)))
    vcov_theta <- tryCatch(solve(opt$hessian), error = function(e) NULL)

  # ---- Standard errors (reporting scale) ------------------------------------
  if (!is.null(vcov_theta)) {
    se_th   <- sqrt(pmax(0, diag(vcov_theta)))
    std_err <- c(
      se_th[idx_e0], se_th[idx_emax0],
      ec50_est * se_th[idx_log_ec50],           # delta method for ec50
      if (p_lin  > 0) se_th[idx_lin]      else NULL,
      if (p_emax > 0) se_th[idx_emax_cov] else NULL,
      if (p_ec50 > 0) se_th[idx_ec50_cov] else NULL
    )
  } else {
    std_err <- rep(NA_real_, n_par)
  }

  # ---- Coefficients (reporting scale) ----------------------------------------
  coef_names <- c(
    "e0", "emax", "ec50",
    if (p_lin  > 0) paste0("linear.",   colnames(Z_lin))  else character(0),
    if (p_emax > 0) paste0("emax.",     colnames(Z_emax)) else character(0),
    if (p_ec50 > 0) paste0("log_ec50.", colnames(Z_ec50)) else character(0)
  )
  coefs <- setNames(
    c(
      th[idx_e0], th[idx_emax0], ec50_est,
      if (p_lin  > 0) th[idx_lin]      else NULL,
      if (p_emax > 0) th[idx_emax_cov] else NULL,
      if (p_ec50 > 0) th[idx_ec50_cov] else NULL
    ),
    coef_names
  )
  names(std_err) <- coef_names

  # ---- Fitted values ---------------------------------------------------------
  emax_i <- th[idx_emax0] +
    if (p_emax > 0) drop(Z_emax %*% th[idx_emax_cov]) else 0
  ec50_i <- exp(th[idx_log_ec50] +
    if (p_ec50 > 0) drop(Z_ec50 %*% th[idx_ec50_cov]) else 0)
  eta    <- th[idx_e0] + emax_i * x / (ec50_i + x) +
    if (p_lin > 0) drop(Z_lin %*% th[idx_lin]) else 0
  fitted <- as.vector(stats::plogis(eta))

  structure(
    list(
      coefficients      = coefs,
      std_err           = std_err,
      vcov_theta        = vcov_theta,
      theta_hat         = th,
      fitted.values     = fitted,
      linear.predictors = as.vector(eta),
      residuals         = y - fitted,
      loglik            = -opt$value,
      AIC               = 2 * n_par + 2 * opt$value,
      converged         = opt$convergence == 0,
      n_par             = n_par,
      x                 = x,
      y                 = y,
      Z_lin             = Z_lin,
      Z_emax            = Z_emax,
      Z_ec50            = Z_ec50,
      idx               = idx,
      call              = match.call()
    ),
    class = "logistic_emax"
  )
}


#' Print a `logistic_emax` object
#'
#' Displays parameter estimates and standard errors, grouped by parameter type.
#'
#' @param x A `logistic_emax` object returned by [fit_logistic_emax()].
#' @param digits Number of decimal places. Default `4`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return `x`, invisibly.
#'
#' @seealso [fit_logistic_emax()], [predict.logistic_emax()]
#'
#' @export
print.logistic_emax <- function(x, digits = 4, ...) {
  fmt <- function(nms, strip_prefix = NULL) {
    tab <- cbind(Estimate = x$coefficients[nms], `Std. Error` = x$std_err[nms])
    if (!is.null(strip_prefix))
      rownames(tab) <- sub(paste0("^", strip_prefix, "\\."), "", nms)
    round(tab, digits)
  }

  cat("Logistic Emax Model\n\n")
  cat("Core parameters:\n")
  print(fmt(c("e0", "emax", "ec50")), ...)

  lin_nm  <- grep("^linear\\.",   names(x$coefficients), value = TRUE)
  emax_nm <- grep("^emax\\.",     names(x$coefficients), value = TRUE)
  ec50_nm <- grep("^log_ec50\\.", names(x$coefficients), value = TRUE)

  if (length(lin_nm) > 0) {
    cat("\nLinear covariates (additive on link scale):\n")
    print(fmt(lin_nm,  "linear"), ...)
  }
  if (length(emax_nm) > 0) {
    cat("\nEmax covariates (additive effect on Emax):\n")
    print(fmt(emax_nm, "emax"), ...)
  }
  if (length(ec50_nm) > 0) {
    cat("\nEC50 covariates (additive effect on log EC50):\n")
    print(fmt(ec50_nm, "log_ec50"), ...)
  }

  cat(
    "\nLog-likelihood:", round(x$loglik, 2),
    "  AIC:", round(x$AIC, 2),
    "  Converged:", x$converged, "\n"
  )
  invisible(x)
}


#' Predict from a `logistic_emax` object
#'
#' Generates predictions on the link, response, or confidence-interval scale.
#' For models with covariates, new covariate matrices can be supplied via
#' `newlinear`, `newemax`, and `newec50`; when these are `NULL` and a new `x`
#' vector is given, covariates are fixed at their training-data column means.
#'
#' Confidence intervals (`type = "ci"`) use the delta method on the link scale,
#' propagating uncertainty through all model parameters including covariate
#' effects.
#'
#' @param object A `logistic_emax` object.
#' @param newdata Optional numeric vector of x values at which to predict.
#'   `NULL` (default) predicts at the original observed values.
#' @param newlinear Optional matrix of linear covariates at the prediction
#'   points. Must have `length(newdata)` rows and the same number of columns as
#'   the training `linear_covs`. When `NULL` and `newdata` is supplied, training
#'   column means are used (with a message).
#' @param newemax Optional matrix of Emax covariates at the prediction points.
#'   Same convention as `newlinear`.
#' @param newec50 Optional matrix of EC50 covariates at the prediction points.
#'   Same convention as `newlinear`.
#' @param type One of `"link"`, `"response"`, or `"ci"`. Default `"link"`.
#' @param level Confidence level for `type = "ci"`. Default `0.95`.
#' @param ... Currently unused; included for S3 compatibility.
#'
#' @return For `type = "link"` or `"response"`, a numeric vector. For
#'   `type = "ci"`, a data frame with columns `x`, `fit`, `lwr`, `upr`.
#'
#' @seealso [fit_logistic_emax()]
#'
#' @examples
#' set.seed(1)
#' dose <- sort(rep(c(0, 5, 10, 25, 50), length.out = 150))
#' eta  <- qlogis(0.10) + 2 * dose / (10 + dose)
#' y    <- rbinom(150, 1, plogis(eta))
#' fit  <- fit_logistic_emax(y, dose)
#'
#' grid <- seq(0, 50, length.out = 100)
#' head(predict(fit, newdata = grid, type = "ci"))
#'
#' @export
predict.logistic_emax <- function(
  object,
  newdata   = NULL,
  newlinear = NULL,
  newemax   = NULL,
  newec50   = NULL,
  type      = c("link", "response", "ci"),
  level     = 0.95,
  ...
) {
  type             <- match.arg(type)
  newdata_provided <- !is.null(newdata)
  x_new            <- if (newdata_provided) as.numeric(newdata) else object$x
  n_new            <- length(x_new)

  # ---- Resolve covariate matrices at prediction points -----------------------
  resolve_cov <- function(new_z, train_z, label) {
    p <- ncol(train_z)
    if (p == 0L) return(matrix(0, nrow = n_new, ncol = 0L))
    if (!is.null(new_z)) {
      m <- as.matrix(new_z)
      if (nrow(m) != n_new)
        stop(label, " must have ", n_new, " rows.")
      if (ncol(m) != p)
        stop(label, " must have ", p, " column(s).")
      return(m)
    }
    if (newdata_provided) {
      message(label, " not supplied; fixing at training column means.")
      return(matrix(colMeans(train_z), nrow = n_new, ncol = p, byrow = TRUE))
    }
    train_z
  }

  Z_lin_new  <- resolve_cov(newlinear, object$Z_lin,  "newlinear")
  Z_emax_new <- resolve_cov(newemax,   object$Z_emax, "newemax")
  Z_ec50_new <- resolve_cov(newec50,   object$Z_ec50, "newec50")

  p_lin  <- ncol(Z_lin_new)
  p_emax <- ncol(Z_emax_new)
  p_ec50 <- ncol(Z_ec50_new)

  # ---- Linear predictor ------------------------------------------------------
  idx <- object$idx
  th  <- object$theta_hat

  emax_i <- th[idx$emax0] +
    if (p_emax > 0) drop(Z_emax_new %*% th[idx$emax_cov]) else 0
  ec50_i <- exp(th[idx$log_ec50] +
    if (p_ec50 > 0) drop(Z_ec50_new %*% th[idx$ec50_cov]) else 0)
  eta    <- th[idx$e0] + emax_i * x_new / (ec50_i + x_new) +
    if (p_lin > 0) drop(Z_lin_new %*% th[idx$lin]) else 0

  if (type == "link")     return(as.vector(eta))
  if (type == "response") return(as.vector(stats::plogis(eta)))

  # ---- CI via delta method ---------------------------------------------------
  V <- object$vcov_theta
  if (is.null(V)) stop("vcov_theta not available; refit with hessian = TRUE.")

  n_par <- object$n_par
  G     <- matrix(0, nrow = n_new, ncol = n_par)
  denom <- ec50_i + x_new

  G[, idx$e0]       <- 1
  G[, idx$emax0]    <- x_new / denom
  G[, idx$log_ec50] <- -emax_i * x_new * ec50_i / denom^2

  # sweep(M, 1, v, "*") multiplies row i of M by v[i]
  if (p_lin  > 0) G[, idx$lin]      <- Z_lin_new
  if (p_emax > 0) G[, idx$emax_cov] <- sweep(Z_emax_new, 1, x_new / denom, `*`)
  if (p_ec50 > 0) G[, idx$ec50_cov] <- sweep(
    Z_ec50_new, 1, -emax_i * x_new * ec50_i / denom^2, `*`
  )

  var_eta <- rowSums((G %*% V) * G)
  se_eta  <- sqrt(pmax(0, var_eta))
  z       <- stats::qnorm(1 - (1 - level) / 2)

  data.frame(
    x   = x_new,
    fit = stats::plogis(eta),
    lwr = stats::plogis(eta - z * se_eta),
    upr = stats::plogis(eta + z * se_eta)
  )
}
