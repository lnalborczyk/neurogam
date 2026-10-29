# Prepare predicted and observed EEG trajectories for plotting

Observations are averaged within each sensor-time location. For a model
contrast A - B, supply \`observed_data\` for A and \`observed_data_b\`
for B. The displayed residual is the difference between observed
\*mean\* and posterior expected value, not a trial-level residual.

## Usage

``` r
prepare_spatiotemporal_plot_data(
  results,
  observed_data,
  observed_data_b = NULL,
  outcome_id = "voltage",
  sensor_id = "channel",
  time_id = "time",
  spatial_id = c("xproj", "yproj"),
  participant_id = NULL
)
```

## Arguments

- results:

  Result from \[analyse_bgam_spatiotemporal()\].

- observed_data:

  Observed data in long format (condition A if relevant).

- observed_data_b:

  Optional observed condition B for a model contrast.

- outcome_id:

  Outcome column of observed data.

- sensor_id:

  Sensor column name in both data and \`results\$predictions\`.

- time_id:

  Numeric time column.

- spatial_id:

  Two projected coordinate column names.

- participant_id:

  Optional participant column for an equal-participant- weighted grand
  average (average within participant first).

## Value

A data frame with standardized sensor, time, x, y, observed, predicted,
lower, upper, residual, cluster_id and cluster_sign columns.
