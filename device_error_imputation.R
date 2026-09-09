## =============================================================================
##  Device-error multiple imputation and Rubin pooling
##
##  Reference implementation accompanying:
##    Casiraghi et al., "Modern Times: Longitudinal Study of Toba/Qom
##    Communities Reveals Delay and Shortening of Sleep in Real-Time",
##    Current Biology.
##
##  This file contains the procedure used to propagate the measurement
##  difference between two actigraphy systems (Philips Respironics Actiwatch
##  with Actiware scoring, and Axivity AX3 with GGIR scoring) into every
##  statistical estimate reported in the paper.
##
##  The file has no dependencies beyond lme4 and, optionally, emmeans, and is
##  free of any project-specific paths. It can be sourced directly:
##
##      source("device_error_imputation.R")
##
##  A runnable worked example using synthetic data is provided separately in
##  demo_synthetic.R, which reproduces the behaviour of the procedure without
##  requiring access to the study data.
##
##  Licence: MIT.
## =============================================================================


## -----------------------------------------------------------------------------
##  1. THE PROBLEM
## -----------------------------------------------------------------------------
##
##  Sleep episodes recorded before 2023 were measured with an Actiwatch; those
##  from 2023 onwards with an AX3. The two systems do not return identical
##  estimates for the same night. In a concurrent-wear validation study, 22
##  volunteers wore both devices simultaneously, yielding 206 paired sleep
##  episodes and hence 206 paired differences per sleep variable.
##
##  Those differences are neither negligible nor normally distributed: they are
##  skewed and heavy-tailed, with mean and median far apart (for sleep duration,
##  mean -17.5 min but median +1.2 min). Two consequences follow:
##
##    (a) Subtracting a single constant offset would be wrong, because it would
##        assume the discrepancy is a fixed bias and would silently treat the
##        corrected values as if they were measured without error.
##
##    (b) Assuming a normal error distribution would be wrong, because the
##        empirical differences are visibly non-normal.
##
##  The procedure below therefore treats the device difference as a source of
##  *missing information* rather than as a fixed correction: we do not know what
##  the pre-2023 records would have read on an AX3, so we impute that quantity
##  repeatedly, drawing from the observed distribution of differences, and
##  propagate the resulting variability into the standard errors.
##
##  The consequence worth stating plainly is that this procedure can only widen
##  confidence intervals relative to an uncorrected analysis. It is conservative
##  by construction: it cannot manufacture a significant result, and every
##  effect reported in the paper survives it.


## -----------------------------------------------------------------------------
##  2. THE PROCEDURE
## -----------------------------------------------------------------------------
##
##  For each of M = 100 imputations m = 1, ..., M:
##
##    Step 1  Draw n values with replacement from the 206 observed paired
##            differences, where n is the number of pre-2023 (Actiwatch)
##            observations in the dataset. Sampling with replacement from the
##            empirical distribution preserves its skew and its tails exactly,
##            with no parametric assumption.
##
##    Step 2  Add one drawn difference to each pre-2023 observation of the
##            outcome, placing the whole series on the AX3 measurement scale.
##            Differences are oriented as (AX3 - Actiwatch), so addition is the
##            correct direction: it maps Actiwatch readings onto the AX3 scale,
##            which is the scale of the most recent data.
##
##    Step 3  Refit the model of interest to this imputed dataset, and compute
##            any marginal means or contrasts required.
##
##  The M sets of estimates are then combined with Rubin's rules (Rubin 1987):
##
##      Qbar = (1/M) * sum_m Q_m                       pooled point estimate
##      Ubar = (1/M) * sum_m U_m                       within-imputation variance
##      B    = (1/(M-1)) * sum_m (Q_m - Qbar)^2        between-imputation variance
##      T    = Ubar + (1 + 1/M) * B                    total variance
##      SE   = sqrt(T)
##      lambda = (1 + 1/M) * B / T                     fraction of missing information
##      nu   = (M - 1) / lambda^2                      degrees of freedom
##
##  Here U_m is the squared standard error of the estimate in imputation m.
##  The term (1 + 1/M) * B is the entire contribution of the device change: it
##  is the extra uncertainty that would vanish if the two devices agreed
##  exactly. Wald statistics and confidence intervals are formed on nu degrees
##  of freedom.
##
##  Note on the degrees-of-freedom formula: nu = (M - 1) / lambda^2 is the
##  original Rubin (1987) expression. It is used here rather than the
##  Barnard-Rubin (1999) small-sample adjustment because the complete-data
##  residual degrees of freedom in this study are very large (thousands of sleep
##  episodes), the regime in which the two expressions coincide. When lambda is
##  small -- that is, when the device change contributes little relative to
##  sampling variability -- nu becomes very large and the reference distribution
##  approaches the normal, which is the intended behaviour.


## -----------------------------------------------------------------------------
##  3. CORE FUNCTIONS
## -----------------------------------------------------------------------------

#' Pool a vector of estimates and their variances using Rubin's rules.
#'
#' This is the single computational kernel of the method. Every pooled quantity
#' reported in the paper -- fixed-effect coefficients, estimated marginal means,
#' and contrasts between years -- is produced by this function.
#'
#' @param estimates Numeric vector of length M: the point estimate from each
#'   imputation.
#' @param variances Numeric vector of length M: the squared standard error of
#'   that estimate in each imputation.
#' @return A one-row data frame with the pooled estimate, standard error,
#'   confidence limits, degrees of freedom, the fraction of missing information,
#'   and the Wald t statistic and p value.
pool_rubin <- function(estimates, variances, conf_level = 0.95) {
  stopifnot(length(estimates) == length(variances), length(estimates) >= 2)
  if (anyNA(estimates) || anyNA(variances)) {
    stop("pool_rubin(): NA in estimates or variances; ",
         "check that every imputation produced a usable fit.")
  }

  M <- length(estimates)

  Q_bar <- mean(estimates)          # pooled point estimate
  U_bar <- mean(variances)          # average within-imputation variance
  B     <- stats::var(estimates)    # between-imputation variance (divides by M-1)
  T_var <- U_bar + (1 + 1 / M) * B  # total variance
  SE    <- sqrt(T_var)

  # Fraction of missing information attributable to the device change. If the
  # two devices agreed exactly, B would be 0, lambda would be 0, and the result
  # would reduce to an ordinary single-fit analysis.
  lambda <- (1 + 1 / M) * B / T_var

  # Guard the degenerate case B = 0 exactly (lambda = 0 gives nu = Inf, which is
  # the correct limit but produces a warning in some pt() implementations).
  df <- if (lambda <= .Machine$double.eps) Inf else (M - 1) / lambda^2

  t_stat <- Q_bar / SE
  p_val  <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
  crit   <- stats::qt(1 - (1 - conf_level) / 2, df = df)

  data.frame(
    estimate  = Q_bar,
    SE        = SE,
    conf.low  = Q_bar - crit * SE,
    conf.high = Q_bar + crit * SE,
    df        = df,
    fmi       = lambda,
    statistic = t_stat,
    p.value   = p_val,
    row.names = NULL
  )
}


#' Pool the fixed effects of a list of fitted lmer models.
#'
#' @param model_list List of M models fitted by lme4::lmer() to the M imputed
#'   datasets. All models must share the same fixed-effect structure.
#' @return A data frame with one row per fixed-effect term.
pool_lmer_fixef <- function(model_list, conf_level = 0.95) {
  stopifnot(length(model_list) >= 2)

  # Matrix of coefficients: rows = terms, columns = imputations.
  coefs <- sapply(model_list, lme4::fixef)
  # Corresponding sampling variances, taken from the diagonal of vcov.
  vars  <- sapply(model_list, function(mod) diag(as.matrix(stats::vcov(mod))))

  if (is.null(dim(coefs))) {
    stop("pool_lmer_fixef(): models appear to have a single fixed effect or ",
         "inconsistent structures across imputations.")
  }
  if (!identical(rownames(coefs), rownames(vars))) {
    stop("pool_lmer_fixef(): coefficient and variance names do not match.")
  }

  out <- do.call(rbind, lapply(seq_len(nrow(coefs)), function(i) {
    pool_rubin(coefs[i, ], vars[i, ], conf_level = conf_level)
  }))
  cbind(term = rownames(coefs), out)
}


#' Pool a set of estimates that are grouped, such as marginal means or contrasts.
#'
#' Used for the emmeans output, where the same quantity is estimated once per
#' imputation within each combination of grouping variables.
#'
#' @param df Data frame stacking the results of all M imputations.
#' @param by Character vector of column names identifying a unique quantity.
#' @param estimate,variance Names of the columns holding the per-imputation
#'   estimate and its squared standard error.
#' @return One pooled row per combination of the `by` columns.
pool_grouped <- function(df, by, estimate = "estimate", variance = "var",
                         conf_level = 0.95) {
  stopifnot(all(c(by, estimate, variance) %in% names(df)))

  key    <- interaction(df[by], drop = TRUE, lex.order = TRUE)
  splits <- split(df, key)

  out <- do.call(rbind, lapply(splits, function(part) {
    pooled <- pool_rubin(part[[estimate]], part[[variance]], conf_level)
    labels <- unique(as.data.frame(part)[, by, drop = FALSE])
    cbind(labels, pooled, row.names = NULL)
  }))
  rownames(out) <- NULL
  out
}


#' Draw one device-error imputation of a dataset.
#'
#' Adds an error drawn from the empirical distribution of paired device
#' differences to every observation measured with the older device.
#'
#' @param data Data frame containing the outcome.
#' @param outcome Name of the outcome column to be shifted.
#' @param deltas Numeric vector of observed paired differences, oriented as
#'   (new device - old device).
#' @param rows_old Integer or logical index of the rows measured with the OLD
#'   device. These are the rows that get shifted onto the new device's scale.
#'   Note the direction carefully: it is the older records that are adjusted,
#'   because the newer device defines the reference scale.
#' @return A copy of `data` with the outcome shifted on the selected rows.
impute_device_error <- function(data, outcome, deltas, rows_old) {
  stopifnot(outcome %in% names(data), length(deltas) > 0)

  if (is.logical(rows_old)) rows_old <- which(rows_old)
  if (length(rows_old) == 0L) {
    warning("impute_device_error(): no rows selected as old-device records; ",
            "the data will be returned unchanged.")
    return(data)
  }

  drawn <- sample(deltas, size = length(rows_old), replace = TRUE)
  data[[outcome]][rows_old] <- data[[outcome]][rows_old] + drawn
  data
}


#' Run a complete multiple-imputation analysis.
#'
#' Convenience wrapper tying the three steps together: impute, refit, pool.
#'
#' @param data Complete dataset.
#' @param formula Model formula passed to lme4::lmer().
#' @param outcome Name of the outcome column (must match the formula's LHS).
#' @param deltas Empirical paired differences (new device - old device).
#' @param rows_old Index of rows measured with the old device.
#' @param M Number of imputations. The paper uses 100.
#' @param seed Random seed. Fixing it makes the analysis exactly reproducible.
#' @param extract Optional function(model, data) returning a data frame of
#'   additional per-imputation quantities to be pooled, for example emmeans
#'   contrasts. It must return the same rows and columns for every imputation,
#'   including a column of estimates and a column of their variances.
#' @return A list with the pooled fixed effects, the pooled extras (if any), the
#'   list of fitted models, and a record of the settings used.
run_mi_analysis <- function(data, formula, outcome, deltas, rows_old,
                            M = 100, seed = 123, extract = NULL,
                            conf_level = 0.95, verbose = TRUE) {

  set.seed(seed)

  models <- vector("list", M)
  extras <- vector("list", M)
  n_singular <- 0L
  n_warned   <- 0L

  for (m in seq_len(M)) {
    imputed <- impute_device_error(data, outcome, deltas, rows_old)

    # Convergence warnings are captured rather than suppressed, so that they can
    # be counted and reported honestly alongside the results.
    fit <- withCallingHandlers(
      lme4::lmer(formula, data = imputed),
      warning = function(w) {
        n_warned <<- n_warned + 1L
        invokeRestart("muffleWarning")
      }
    )

    if (lme4::isSingular(fit)) n_singular <- n_singular + 1L

    models[[m]] <- fit
    if (!is.null(extract)) extras[[m]] <- extract(fit, imputed)

    if (verbose && m %% 10 == 0) message("  imputation ", m, " / ", M)
  }

  pooled_fixef <- pool_lmer_fixef(models, conf_level = conf_level)

  pooled_extra <- NULL
  if (!is.null(extract)) pooled_extra <- do.call(rbind, extras)

  if (verbose) {
    message("Completed ", M, " imputations.")
    message("  singular fits: ", n_singular, " / ", M)
    message("  fits emitting warnings: ", n_warned, " / ", M)
  }

  list(
    fixed_effects = pooled_fixef,
    per_imputation = pooled_extra,
    models = models,
    settings = list(M = M, seed = seed, outcome = outcome,
                    n_old_rows = length(if (is.logical(rows_old))
                      which(rows_old) else rows_old),
                    n_deltas = length(deltas),
                    n_singular = n_singular, n_warned = n_warned)
  )
}


## -----------------------------------------------------------------------------
##  4. REFERENCES
## -----------------------------------------------------------------------------
##
##  Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
##  Barnard, J. and Rubin, D.B. (1999). Small-sample degrees of freedom with
##      multiple imputation. Biometrika 86, 948-955.
##  Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting linear
##      mixed-effects models using lme4. J. Stat. Softw. 67, 1-48.
##  Lenth, R.V. emmeans: Estimated Marginal Means, aka Least-Squares Means.
