## =============================================================================
##  Worked example: device-error multiple imputation using the REAL measured
##  device differences
##
##  This is the companion to demo_synthetic.R. The difference is where the
##  device error comes from:
##
##    demo_synthetic.R   simulated differences drawn from a skewed mixture
##    demo_real.R        the 206 measured differences in data/device_comparison.csv
##
##  The study cohort is still synthetic, because the participant data is under
##  controlled access. But the measurement error applied to it is exactly the
##  error that was measured between the two devices and used in the paper, so
##  this script shows the real magnitude and shape of the correction.
##
##  Run from this directory:  Rscript demo_real.R
##
##  Requires: lme4.
## =============================================================================

source("device_error_imputation.R")

suppressPackageStartupMessages(library(lme4))

set.seed(2024)


## -----------------------------------------------------------------------------
##  1. The measured device differences
## -----------------------------------------------------------------------------

dev <- read.csv("data/device_comparison.csv", stringsAsFactors = FALSE)

stopifnot(nrow(dev) == 206L, length(unique(dev$volunteer_id)) == 22L)

deltas_onset    <- dev$diff_onset       # hours, AX3 - Actiwatch
deltas_duration <- dev$diff_duration

describe <- function(label, v) {
  cat(sprintf("  %-16s n=%3d  mean %+6.1f  median %+6.1f  SD %5.1f  ",
              label, length(v), mean(v) * 60, median(v) * 60, sd(v) * 60))
  cat(sprintf("2.5%%..97.5%% [%+.1f, %+.1f] min\n",
              quantile(v, .025) * 60, quantile(v, .975) * 60))
}

cat("Measured differences (AX3 - Actiwatch), minutes:\n")
describe("sleep onset",    deltas_onset)
describe("sleep duration", deltas_duration)

# The distance between mean and median is the whole reason for resampling.
cat(sprintf("\n  Mean-minus-median, onset:    %+.1f min\n",
            (mean(deltas_onset) - median(deltas_onset)) * 60))
cat(sprintf("  Mean-minus-median, duration: %+.1f min\n",
            (mean(deltas_duration) - median(deltas_duration)) * 60))
cat("  A constant-offset correction would apply the mean to every night,\n")
cat("  even though most nights sit near the median.\n")


## -----------------------------------------------------------------------------
##  2. A synthetic cohort with the study's structure
## -----------------------------------------------------------------------------
##  120 participants, two communities, campaigns from 2012 to 2024, ~15 nights
##  each per campaign, and a built-in true delay of +0.09 h/year in sleep onset.

n_part    <- 120
campaigns <- c(2012, 2013, 2014, 2016, 2017, 2018, 2023, 2024)
true_slope <- 0.09

participants <- data.frame(
  uid       = sprintf("P%03d", seq_len(n_part)),
  community = rep(c("rural", "urban"), length.out = n_part),
  intercept = rnorm(n_part, 0, 0.6)
)

rows <- do.call(rbind, lapply(seq_len(n_part), function(i) {
  # Unbalanced follow-up: each participant appears in a random subset.
  mine   <- sort(sample(campaigns, size = sample(2:5, 1)))
  nights <- 15
  data.frame(
    uid       = participants$uid[i],
    community = participants$community[i],
    year      = rep(mine, each = nights),
    b_i       = participants$intercept[i]
  )
}))

rows$onset_true <- 22.0 +
  true_slope * (rows$year - 2012) +
  ifelse(rows$community == "urban", 0.40, 0) +
  rows$b_i +
  rnorm(nrow(rows), 0, 1.0)


## -----------------------------------------------------------------------------
##  3. What was actually observed
## -----------------------------------------------------------------------------
##  Records before 2023 came off the Actiwatch. To simulate that, subtract a
##  real measured difference from each: if the AX3 reads later than the
##  Actiwatch by d, then the Actiwatch would have recorded truth - d.

rows$old_device <- rows$year < 2023

rows$onset_obs <- rows$onset_true
rows$onset_obs[rows$old_device] <-
  rows$onset_true[rows$old_device] -
  sample(deltas_onset, sum(rows$old_device), replace = TRUE)

cat(sprintf("\nSynthetic cohort: %d participants, %d nights, %d on the old device.\n",
            n_part, nrow(rows), sum(rows$old_device)))


## -----------------------------------------------------------------------------
##  4. Three analyses
## -----------------------------------------------------------------------------

form <- onset_obs ~ year_c + community + (1 | uid)
rows$year_c <- rows$year - 2012

ctrl <- lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))

# (a) Naive: pretend the instrument never changed.
fit_naive <- lmer(form, data = rows, REML = FALSE, control = ctrl)
est_naive <- summary(fit_naive)$coefficients["year_c", ]

# (b) Constant offset: add the mean difference to the old records, then analyse
#     as though the corrected values had been measured exactly.
rows_const <- rows
rows_const$onset_obs[rows_const$old_device] <-
  rows_const$onset_obs[rows_const$old_device] + mean(deltas_onset)
fit_const <- lmer(form, data = rows_const, REML = FALSE, control = ctrl)
est_const <- summary(fit_const)$coefficients["year_c", ]

# (c) Multiple imputation from the measured empirical distribution.
mi <- run_mi_analysis(
  data     = rows,
  formula  = form,
  outcome  = "onset_obs",
  deltas   = deltas_onset,
  rows_old = rows$old_device,
  M        = 100,
  seed     = 123,
  verbose  = FALSE
)
est_mi <- mi$fixed_effects[mi$fixed_effects$term == "year_c", ]
stopifnot(nrow(est_mi) == 1L)


## -----------------------------------------------------------------------------
##  5. Results
## -----------------------------------------------------------------------------

cat(sprintf("\n\nEffect of year on sleep onset (true value %+.4f h/year)\n\n",
            true_slope))
cat(sprintf("  %-26s %9s %8s   %-22s %7s\n",
            "", "estimate", "SE", "95% CI", "width"))

show <- function(label, est, se, lo, hi)
  cat(sprintf("  %-26s %+9.4f %8.4f   [%+.4f, %+.4f] %7.4f\n",
              label, est, se, lo, hi, hi - lo))

show("(a) naive, uncorrected", est_naive[1], est_naive[2],
     est_naive[1] - 1.96 * est_naive[2], est_naive[1] + 1.96 * est_naive[2])
show("(b) constant mean offset", est_const[1], est_const[2],
     est_const[1] - 1.96 * est_const[2], est_const[1] + 1.96 * est_const[2])
show("(c) multiple imputation", est_mi$estimate, est_mi$SE,
     est_mi$conf.low, est_mi$conf.high)

cat(sprintf("\n  Fraction of missing information attributable to the device change: %.1f%%\n",
            100 * est_mi$fmi))
cat(sprintf("  Interval widening from (b) to (c): %.0f%%\n",
            100 * (((est_mi$conf.high - est_mi$conf.low) /
                    (2 * 1.96 * est_const[2])) - 1)))

cat("\n")
cat("Reading the table:\n")
cat("  (a) is biased. The device difference is absorbed into the year term,\n")
cat("      because the change of instrument coincides with the passage of time.\n")
cat("  (b) recovers the point estimate but reports the same SE as (a). It\n")
cat("      treats the corrected values as if measured without error.\n")
cat("  (c) recovers the point estimate AND widens the interval to reflect what\n")
cat("      is actually known about the instrument. This is the conservative\n")
cat("      option: it can only make an effect harder to detect, never easier.\n")

cat("\nDone.\n")
