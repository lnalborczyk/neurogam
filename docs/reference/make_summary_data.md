# Summarise trial-level data for `testing_through_time()`

Summarise trial-level data for
[`testing_through_time()`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md)

## Usage

``` r
make_summary_data(
  data,
  participant_id = "participant",
  outcome_id = "eeg",
  outcome_sd = NULL,
  time_id = "time",
  predictor_id = NA,
  trials_id = NULL,
  family = gaussian(),
  multilevel = c("summary", "group"),
  na_rm = TRUE,
  continuous_agg = c("mean", "first")
)
```
