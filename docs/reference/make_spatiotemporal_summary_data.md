# Create participant-level EEG summaries for spatial-temporal fitting

Create participant-level EEG summaries for spatial-temporal fitting

## Usage

``` r
make_spatiotemporal_summary_data(
  data,
  participant_id,
  outcome_id,
  outcome_sem = NULL,
  time_id = "time",
  sensor_id = "channel",
  space_id = c("xproj", "yproj"),
  predictor_id = NA_character_,
  use_se = TRUE
)
```

## Arguments

- data:

  Trial-level or previously summarized data frame.

- participant_id, :

  outcome_id, outcome_sem, time_id, sensor_id, space_id, predictor_id
  Source column names; \`predictor_id\` may be \`NA\`.

- use_se:

  Whether the model will use observation-level SEM.

## Value

Data frame with standardized model columns.
