# Device-error multiple imputation and Rubin pooling

Reference implementation of the measurement-error correction used in:

> Casiraghi et al. *Modern Times: Longitudinal Study of Toba/Qom Communities
> Reveals Delay and Shortening of Sleep in Real-Time.* Current Biology, 2026.

## What this is

Sleep episodes in this study were recorded with two different actigraphy
systems: a Philips Respironics Actiwatch scored with Actiware 6 (2012–2018) and
an Axivity AX3 scored with GGIR 3.2.6 (2023–2024). Because the change of instrument
coincides with the passage of time, any difference between the two systems is
perfectly confounded with the temporal trend that the study set out to measure.

This code implements the procedure used to prevent that confound from being
mistaken for a real effect. It treats the device difference as *missing
information* — we do not know what the earlier records would have read on an
AX3 — and imputes it repeatedly from the empirical distribution of paired
differences measured in a concurrent-wear validation study, then combines the
results with Rubin's rules so that the device uncertainty is carried through
into every standard error and confidence interval.

The important property is that **the procedure is conservative by
construction**. Adding between-imputation variance can only widen intervals
relative to an uncorrected analysis. It cannot manufacture a significant
result. Every effect reported in the paper survives it.

## Files

| File | Contents |
|---|---|
| `device_error_imputation.R` | The method. Documented functions, no project-specific paths, no external data. |
| `demo_synthetic.R` | Runnable worked example on entirely synthetic data. Requires no external files. |
| `demo_real.R` | The same worked example driven by the **measured** device differences in `data/`. |
| `data/device_comparison.csv` | The 206 paired concurrent-wear differences used in the paper. See `data/README.md`. |

## Licensing

The code is under the MIT Licence (`LICENSE`). The dataset in `data/` is under
**CC BY 4.0**; see `data/README.md`.

## Running the demos

```r
# from this directory
Rscript demo_synthetic.R    # synthetic cohort, synthetic device error
Rscript demo_real.R         # synthetic cohort, real measured device error
```

Requires `lme4`; `emmeans` is optional and used only for the final section.

The demo generates a synthetic dataset with the study's structure — two
communities, unbalanced longitudinal follow-up, a device switch partway
through, and a known true effect of +0.090 h/year — and compares three
analyses. Representative output:

```
                             estimate  SE       95% CI                 width
  (a) naive, uncorrected     +0.0787   0.0083   [+0.0624, +0.0950]   0.0326
  (b) constant mean offset   +0.0892   0.0083   [+0.0729, +0.1055]   0.0326
  (c) multiple imputation    +0.0890   0.0100   [+0.0695, +0.1086]   0.0391
```

This makes the argument concretely:

1. The **naive** analysis is biased — the device difference leaks into the
   year term (0.079 against a true 0.090).
2. A **constant offset** recovers the point estimate but leaves the standard
   error untouched. It treats corrected values as though they had been
   measured exactly, and so understates what is actually known.
3. **Multiple imputation** recovers the point estimate *and* widens the
   interval by about 20% to reflect the device uncertainty.

## The method

For each of M = 100 imputations:

1. Draw *n* values with replacement from the 206 observed paired differences,
   where *n* is the number of Actiwatch-era observations. Sampling from the
   empirical distribution preserves its skew and heavy tails exactly, with no
   parametric assumption. This matters: for sleep duration the observed
   differences have a mean of −17.5 min but a median of +1.2 min, so a normal
   approximation would misrepresent them.
2. Add one drawn difference to each Actiwatch-era observation, placing the
   whole series on the AX3 scale. Differences are oriented as
   (AX3 − Actiwatch), so the *older* records are the ones adjusted — the newer
   device defines the reference scale.
3. Refit the mixed model and compute any marginal means or contrasts required.

The M results are then combined by Rubin's rules (Rubin 1987):

```
Qbar   = mean of the M estimates                    pooled estimate
Ubar   = mean of the M sampling variances           within-imputation variance
B      = variance of the M estimates                between-imputation variance
T      = Ubar + (1 + 1/M) * B                       total variance
lambda = (1 + 1/M) * B / T                          fraction of missing information
nu     = (M - 1) / lambda^2                         degrees of freedom
```

The term `(1 + 1/M) * B` is the entire contribution of the device change: it is
the extra uncertainty that would disappear if the two systems agreed exactly.
`lambda` is directly interpretable as the proportion of the total uncertainty
in a given estimate that is due to the instrument change, and is reported by
`pool_rubin()` as `fmi`.

The degrees-of-freedom expression is the original Rubin (1987) form rather than
the Barnard–Rubin (1999) small-sample adjustment, because the complete-data
residual degrees of freedom here are in the thousands of sleep episodes — the
regime in which the two coincide.

## Functions

- **`pool_rubin(estimates, variances)`** — the computational kernel. Every
  pooled quantity in the paper passes through this one function.
- **`pool_lmer_fixef(model_list)`** — pools the fixed effects of M fitted
  `lmer` models.
- **`pool_grouped(df, by, ...)`** — pools grouped quantities such as estimated
  marginal means or contrasts.
- **`impute_device_error(data, outcome, deltas, rows_old)`** — one draw.
- **`run_mi_analysis(...)`** — impute, refit and pool in one call. Counts
  singular fits and convergence warnings rather than suppressing them.

## Reproducibility

Seeds are fixed (`seed = 123` in the paper's analyses), so results are exactly
reproducible. `run_mi_analysis()` returns the settings actually used, together
with counts of singular fits and warnings, in its `settings` element.

## Data

The procedure needs only a vector of paired device differences. In this study
those are the 206 concurrent-wear differences from 22 volunteers, **included in
this deposit** as `data/device_comparison.csv`.

| Difference | Column | Mean | Median |
|---|---|---|---|
| Sleep start | `diff_onset` | +16.5 min | +1.4 min |
| Sleep end | `diff_offset` | -0.9 min | +2.6 min |
| Sleep duration | `diff_duration` | -17.5 min | +1.2 min |
| Mid-sleep | `diff_midsleep` | (not used in the paper) | |

All differences are in **hours**, oriented AX3 minus Actiwatch. See
`data/README.md` for the full schema, the two time conventions that matter,
and how the file was produced.

`demo_synthetic.R` requires no data at all; `demo_real.R` reads this file.

The participant data from the Toba/Qom communities is *not* part of this
deposit. It is under controlled access because of the vulnerability of the
communities studied; see the data availability statement in the article.

## References

- Rubin, D.B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
- Barnard, J. & Rubin, D.B. (1999). Small-sample degrees of freedom with
  multiple imputation. *Biometrika* 86, 948–955.
- Bates, D., Mächler, M., Bolker, B. & Walker, S. (2015). Fitting linear
  mixed-effects models using lme4. *J. Stat. Softw.* 67, 1–48.

## Licence

Code: MIT (`LICENSE`).
Data (`data/device_comparison.csv`): CC BY 4.0.
