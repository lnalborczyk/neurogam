# Changelog

## neurogam 0.0.3

### New features

- Implementing the
  [`recommend_k()`](https://lnalborczyk.github.io/neurogam/reference/recommend_k.md)
  function, running model comparison with varying basis dimension (k).
- Implementing a
  [`ppc()`](https://lnalborczyk.github.io/neurogam/reference/ppc.md)
  method for `cluster_results` objects.
- Implementing the
  [`check_residual_autocorrelation()`](https://lnalborczyk.github.io/neurogam/reference/check_residual_autocorrelation.md)
  function to visualise auto-correlation of the residuals.
- New
  [`plot_sensors()`](https://lnalborczyk.github.io/neurogam/reference/plot_sensors.md)
  and
  [`plot_eeg()`](https://lnalborczyk.github.io/neurogam/reference/plot_eeg.md)
  functions to visualise EEG sensors grid, EEG raw data, and GAM
  predictions.
- Adding support for returning participant-level onsets in
  [`testing_through_time()`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
  [`find_clusters()`](https://lnalborczyk.github.io/neurogam/reference/find_clusters.md),
  and related methods such as
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- Adding support for [`binomial()`](https://rdrr.io/r/stats/family.html)
  responses.
- Supporting continuous predictors in
  [`testing_through_time()`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).

### Other changes

- Improved documentation for
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.
- Removed the `multilevel = "full"` option (too slow).
- Now returns clusters with both positive and negative signs in
  [`find_clusters()`](https://lnalborczyk.github.io/neurogam/reference/find_clusters.md).

### Bug fixes

- Fixing aberrant PPCs in the
  [`ppc()`](https://lnalborczyk.github.io/neurogam/reference/ppc.md)
  method (and simplifying arguments with
  `ppc_type = c("group", "participant")`).
- Allowing “negative” clusters in `cluster_results`
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## neurogam 0.0.2

### New features

- Allowing 3 different models to be fitted: full GAMM, GAMM with summary
  statistics (recommended), or group-level GAM.
- Adding further support for presence or absence of predictor (e.g.,
  group, condition). When `predictor_id = NA`, `neurogam` now tests
  signal against 0 through time.
- Implementing [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for
  `cluster_results` objects.
- Improved plotting: now plotting the GAM predictions with raw data and
  improved clusters aesthetics.

### Other changes

- Improved functions documentation and new package website.
- Factoring posterior odds computation within internal functions.

### Bug fixes

- Fixing group-level posterior predictions when multilevel is “full” or
  “summary”.

## neurogam 0.0.1

- Pushing the first version of `neurogam`.
