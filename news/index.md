# Changelog

## neurogam 0.0.3

### New features

- Implementing the
  [`recommend_k()`](https://lnalborczyk.github.io/neurogam/reference/recommend_k.md)
  function, running model comparison with varying basis dimension (k).
- Implementing a
  [`ppc()`](https://lnalborczyk.github.io/neurogam/reference/ppc.md)
  method for `cluster_results` objects.

### Other changes

- Improved documentation for
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.

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
