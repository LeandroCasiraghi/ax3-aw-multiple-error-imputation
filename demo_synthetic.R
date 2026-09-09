## =============================================================================
##  Worked example: device-error multiple imputation on synthetic data
##
##  This script is self-contained. It generates a synthetic dataset with the
##  same structure as the study design and runs the full procedure on it, so
##  that the method can be inspected and executed without access to the
##  participant data.
##
##  Run with:  Rscript demo_synthetic.R
##
##  Requires: lme4. (emmeans is used only for the optional contrast section.)
## =============================================================================

source("device_error_imputation.R")

suppressPackageStartupMessages(library(lme4))

set.seed(2024)


## -----------------------------------------------------------------------------
##  1. Synthetic study data
## -----------------------------------------------------------------------------
##  120 participants in two communities, each measured in several campaigns
##  between 2012 and 2024, with roughly 15 nights per participant per campaign.
##  A true delay in sleep onset of +0.09 h per year is built in, matching the
##  effect size reported in the paper, so that the demonstration is realistic.

n_id      <- 120
campaigns <- c(2012, 2014, 2016, 2018, 2023, 2024)
true_slope <- 0.09   # hours of delay per year -- the quantity of interest

participants <- data.frame(
  uid   = sprintf("P%03d", seq_len(n_id)),
  group = rep(c("Rural", "Urban"), length.out = n_id),
  # participant-level random intercept
  b0    = rnorm(n_id, 0, 0.55)
)

sim <- do.call(rbind, lapply(seq_len(n_id), function(i) {
  # each participant attends a random subset of campaigns (longitudinal but
  # unbalanced, as in the real study)
  attended <- sort(sample(campaigns, size = sample(1:4, 1)))
  do.call(rbind, lapply(attended, function(yr) {
    n_nights <- 15
    data.frame(
      uid   = participants$uid[i],
      group = participants$group[i],
      year  = yr,
      b0    = participants$b0[i],
      night = seq_len(n_nights)
    )
  }))
}))

sim$ym <- sim$year - 2018.5          # centred, as in the paper

# True sleep onset on the AX3 measurement scale.
sim$onset_true <- 23.1 + sim$b0 + true_slope * sim$ym + rnorm(nrow(sim), 0, 0.85)


## -----------------------------------------------------------------------------
##  2. Synthetic device-validation sample
## -----------------------------------------------------------------------------
##  206 paired differences, as in the concurrent-wear validation study. They are
##  deliberately generated as a skewed, heavy-tailed mixture rather than a
##  normal distribution, because that is the salient feature of the real data
##  and the reason a parametric correction would be inadequate.

n_pairs <- 206
deltas <- c(
  rnorm(round(n_pairs * 0.75), mean =  0.02, sd = 0.20),  # bulk: near agreement
  rnorm(round(n_pairs * 0.25), mean = -0.55, sd = 0.90)   # tail: large disagreements
)
deltas <- deltas[seq_len(n_pairs)]

cat("Synthetic validation sample (differences: AX3 - Actiwatch, hours)\n")
cat(sprintf("  n      = %d\n", length(deltas)))
cat(sprintf("  mean   = %+.3f h (%+.1f min)\n", mean(deltas), mean(deltas) * 60))
cat(sprintf("  median = %+.3f h (%+.1f min)\n", median(deltas), median(deltas) * 60))
cat(sprintf("  SD     =  %.3f h\n", sd(deltas)))
cat("  -> mean and median differ substantially: the distribution is skewed,\n")
cat("     which is precisely why a single constant offset is not adequate.\n\n")


## -----------------------------------------------------------------------------
##  3. Construct the observed data
## -----------------------------------------------------------------------------
##  Records before 2023 were taken with the OLD device, so what was actually
##  observed for them is the truth minus a device difference. Records from 2023
##  onwards are already on the new device's scale and are observed directly.

rows_old <- which(sim$year < 2023)

sim$onset <- sim$onset_true
sim$onset[rows_old] <- sim$onset_true[rows_old] -
  sample(deltas, length(rows_old), replace = TRUE)

cat(sprintf("Dataset: %d observations, %d participants, %d pre-2023 records\n\n",
            nrow(sim), length(unique(sim$uid)), length(rows_old)))


## -----------------------------------------------------------------------------
##  4. Three analyses compared
## -----------------------------------------------------------------------------

model_formula <- onset ~ group * ym + (1 | uid)

## extend the limits for pbkrtest and lmerTest to avoid warnings
emmeans::emm_options(pbkrtest.limit = 5000)
emmeans::emm_options(lmerTest.limit = 5000)

## (a) Naive: ignore the device change entirely.
fit_naive <- lmer(model_formula, data = sim)
s_naive   <- summary(fit_naive)$coefficients["ym", ]
ci_naive  <- c(s_naive[1] - 1.96 * s_naive[2], s_naive[1] + 1.96 * s_naive[2])

## (b) Constant offset: shift the old records by the mean difference and then
##     analyse as if the corrected values were measured without error.
sim_offset <- sim
sim_offset$onset[rows_old] <- sim_offset$onset[rows_old] + mean(deltas)
fit_offset <- lmer(model_formula, data = sim_offset)
s_offset   <- summary(fit_offset)$coefficients["ym", ]
ci_offset  <- c(s_offset[1] - 1.96 * s_offset[2], s_offset[1] + 1.96 * s_offset[2])

## (c) Multiple imputation of the device error, pooled by Rubin's rules.
cat("Running multiple imputation (M = 100)...\n")
res_mi <- run_mi_analysis(
  data     = sim,
  formula  = model_formula,
  outcome  = "onset",
  deltas   = deltas,
  rows_old = rows_old,
  M = 100, seed = 123, verbose = FALSE
)
mi_ym <- res_mi$fixed_effects[res_mi$fixed_effects$term == "ym", ]


## -----------------------------------------------------------------------------
##  5. Results
## -----------------------------------------------------------------------------

cat("\n", strrep("=", 72), "\n", sep = "")
cat("Effect of year on sleep onset (hours/year). True value: ",
    sprintf("%+.3f\n", true_slope), sep = "")
cat(strrep("=", 72), "\n\n", sep = "")

row_fmt <- function(label, est, se, lo, hi) {
  cat(sprintf("  %-26s %+.4f   %.4f   [%+.4f, %+.4f]   %.4f\n",
              label, est, se, lo, hi, hi - lo))
}
cat(sprintf("  %-26s %-9s %-8s %-22s %s\n",
            "", "estimate", "SE", "95% CI", "width"))
row_fmt("(a) naive, uncorrected", s_naive[1], s_naive[2], ci_naive[1], ci_naive[2])
row_fmt("(b) constant mean offset", s_offset[1], s_offset[2], ci_offset[1], ci_offset[2])
row_fmt("(c) multiple imputation", mi_ym$estimate, mi_ym$SE,
        mi_ym$conf.low, mi_ym$conf.high)

cat(sprintf("\n  Fraction of missing information attributable to the device change: %.1f%%\n",
            100 * mi_ym$fmi))
cat(sprintf("  Degrees of freedom after pooling: %.0f\n", mi_ym$df))
cat(sprintf("  CI widening, MI relative to constant offset: %+.1f%%\n",
            100 * ((mi_ym$conf.high - mi_ym$conf.low) /
                     (ci_offset[2] - ci_offset[1]) - 1)))

cat("\n", strrep("-", 72), "\n", sep = "")
cat("Points this demonstrates:\n\n")
cat("  1. The naive analysis is biased: the device difference is absorbed into\n")
cat("     the year term, because the change of device coincides with the passage\n")
cat("     of time.\n\n")
cat("  2. The constant-offset correction recovers the point estimate but reports\n")
cat("     a standard error that is too small. It presents corrected values as if\n")
cat("     they had been measured exactly, so its confidence interval is not\n")
cat("     honest about what is actually known.\n\n")
cat("  3. Multiple imputation recovers the point estimate AND widens the interval\n")
cat("     to reflect the device uncertainty. This is the conservative choice: it\n")
cat("     can only ever make an effect harder to detect, never easier.\n")
cat(strrep("-", 72), "\n", sep = "")


## -----------------------------------------------------------------------------
##  6. Optional: pooling grouped quantities such as marginal means
## -----------------------------------------------------------------------------
##  In the paper, estimated marginal means and between-year contrasts are pooled
##  with the same kernel via pool_grouped(). This section shows the pattern; it
##  is skipped if emmeans is unavailable.

if (requireNamespace("emmeans", quietly = TRUE)) {
  cat("\nPooled marginal means by group and year (via pool_grouped):\n\n")

  stacked <- do.call(rbind, lapply(res_mi$models, function(fit) {
    emm <- emmeans::emmeans(fit, ~ ym | group,
                            at = list(ym = c(2014.5, 2024.5) - 2018.5))
    d <- as.data.frame(emm)
    data.frame(group = d$group, ym = d$ym,
               estimate = d$emmean, var = d$SE^2)
  }))

  pooled_emm <- pool_grouped(stacked, by = c("group", "ym"))
  pooled_emm$year <- pooled_emm$ym + 2018.5

  print(data.frame(
    group = pooled_emm$group,
    year  = pooled_emm$year,
    onset = round(pooled_emm$estimate, 2),
    CI    = sprintf("[%.2f, %.2f]", pooled_emm$conf.low, pooled_emm$conf.high)
  ), row.names = FALSE)
} else {
  cat("\n(emmeans not installed; skipping the marginal-means section.)\n")
}

cat("\nDone.\n")
