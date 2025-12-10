# Changelog

## neurogam 0.0.2

### New features

- Improved plotting: now plotting the GAM-smoothed signal with raw data
  and improved clusters aesthetics.
- Allowing 3 different models to be fitted: full GAMM, GAMM with summary
  statistics (recommended), or group-level GAM.
- Adding further support for presence or absence of predictor (e.g.,
  group, condition). When predictor_id=NA, neurogam now tests signal
  against 0 through time.

### Other changes

- Improved functions documentation and package website
- Factoring posterior odds computation within internal functions

### Bug fixes

- Fixing group-level posterior predictions when multilevel is “full” or
  “summary”.

## neurogam 0.0.1

- Pushing the first version of `neurogam`.
