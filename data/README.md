# Device comparison dataset

Paired sleep estimates from 206 nights on which 22 volunteers wore a Philips
Respironics Actiwatch and an Axivity AX3 accelerometer **at the same time**.

These are the empirical paired differences that the multiple-imputation code in
this deposit resamples from. They are released so that the measurement-error
correction reported in the paper can be inspected, reproduced and criticised
independently of the participant data, which is under controlled access.

## File

`device_comparison.csv` — 206 rows, one per volunteer-night.

| Column | Units | Meaning |
|---|---|---|
| `volunteer_id` | — | Pseudonym, `V01`–`V22`. Assigned in order of first recording. |
| `night_date` | date | Calendar date on which the sleep episode began. |
| `onset_actiwatch` | decimal hours | Sleep start, Actiwatch / Actiware. |
| `onset_ax3` | decimal hours | Sleep start, AX3 / GGIR. |
| `offset_actiwatch` | decimal hours | Sleep end, Actiwatch / Actiware. |
| `offset_ax3` | decimal hours | Sleep end, AX3 / GGIR. |
| `duration_actiwatch` | hours | `offset_actiwatch − onset_actiwatch`. |
| `duration_ax3` | hours | `offset_ax3 − onset_ax3`. |
| `diff_onset` | hours | `onset_ax3 − onset_actiwatch`. |
| `diff_offset` | hours | `offset_ax3 − offset_actiwatch`. |
| `diff_duration` | hours | `duration_ax3 − duration_actiwatch`. |
| `diff_midsleep` | hours | Midpoint AX3 − midpoint Actiwatch. |

### Two conventions that will bite you if you miss them

**Times do not wrap at midnight.** Clock times are expressed as decimal hours
measured from 00:00 on `night_date`, so an onset at 00:45 the following morning
is `24.75`, not `0.75`. Wake times are therefore almost always > 24. This is
deliberate: it makes onsets on either side of midnight directly comparable and
makes `duration = offset − onset` valid without a correction term.

**Differences are AX3 minus Actiwatch.** The AX3 defines the reference scale,
because it is the instrument used in the most recent data collection. A
positive `diff_onset` means the AX3 placed sleep onset *later* than the
Actiwatch did. When these values are used to correct a historical series, it is
the older Actiwatch-era records that are shifted onto the AX3 scale.

## How the data were produced

Volunteers wore both devices concurrently on the non-dominant wrist. The two
recordings were then processed through the two complete pipelines exactly as
they were in the main study — no attempt was made to harmonise them — because
the quantity of interest is the difference between *systems*, not between
sensors:

- **Actiwatch** → Actiware, using the manufacturer's rest-interval detection.
- **AX3** → GGIR v3.2.6 for R. Devices were set to 50 Hz with a range of 8 g.
  Raw acceleration was summarised with the default Euclidean Norm Minus One
  (ENMO) metric, and sleep was detected with the implemented van Hees
  algorithm, which identifies sustained immobility bouts from changes in arm
  angle. A minimum of 16 hours of wear time was required for each day.

The AX3 side of this comparison used **the same device model, the same
acquisition settings and the same GGIR version (3.2.6)** as the 2023–2024
participant data reported in the article. The comparison therefore measures
the difference between the two *systems* as actually deployed, not between two
differently configured analyses.

Records were then matched on volunteer and calendar date, and filtered exactly
as in the main analysis:

1. Mean sleep onset (of the two devices) between 19:30 and 05:30.
2. Mean sleep duration within 2.5 median absolute deviations of the median.

The 206 rows here are what survived. Recordings span 2023-03-17 to 2025-02-20;
volunteers contributed between 2 and 14 nights each (median 11). Because the
volunteers were recruited over that period rather than measured in a single
block, the dates are spread out; device model, acquisition settings and
scoring software version were unchanged throughout.

## Summary of the differences

Values in minutes, AX3 − Actiwatch:

| Quantity | Mean | Median |
|---|---|---|
| Sleep start | +16.5 | +1.4 |
| Sleep end | −0.9 | +2.6 |
| Sleep duration | −17.5 | +1.2 |

The gap between the means and the medians is the important feature of this
dataset, and the reason the correction is done by resampling rather than by
subtracting a constant. The distributions are strongly skewed: on most nights
the two systems agree closely, but a minority of nights — typically ones where
the AX3 pipeline placed sleep onset much later — pull the mean well away from
the centre of the distribution. A correction based on a mean offset, or on any
normal approximation, would misrepresent this. Resampling from the observed
values reproduces the skew and the tails exactly, at no distributional cost.

Note also that the mean difference in sleep duration (−17.5 min) and the mean
difference in sleep onset (+16.5 min) are of comparable size and opposite sign,
while sleep end barely differs. The systems disagree mainly about **when sleep
begins**, and the duration difference follows almost entirely from that.

## Participants

Twenty-two adult volunteers, none of them members of the Toba/Qom communities
studied in the paper. They took part in a separate methodological study of
actigraphic comparison approved by the Ethics Committee of Universidad Nacional
de Quilmes (CE-UNQ Nº 6/2023), and all gave informed consent. **The open
release of these de-identified data is covered by that approval.**

Volunteer demographics (age, sex) are reported in the article and are not
included here: with only 22 individuals, the combination of age, sex and a
dated multi-night recording schedule would be more identifying than the
analysis requires. Nothing in the measurement-error correction uses them.

Identifiers are pseudonyms. No name, birthdate, device serial number or
location is present in this file.

## Licence

`device_comparison.csv` is released under **Creative Commons Attribution 4.0
International (CC BY 4.0)**. The code in the parent directory is under the MIT
Licence; see `../LICENSE`.

If you use this dataset, please cite the article — see `../CITATION.cff`.

## Provenance

Generated by `make_device_comparison_data.R` from the study's internal working
file. That script pseudonymises identifiers, verifies that the differences are
internally consistent, runs a disclosure audit, and re-derives the six summary
statistics quoted in the article.

The builder is deliberately **not** included in this deposit: it reads the
internal file and names the volunteers, so shipping it would defeat the
pseudonymisation it performs. It is retained in the authors' repository and is
available on request.
