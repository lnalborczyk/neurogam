# Construct metadata for a neurogam result

Records the analysis specification that is known when a result is
created. This function does not infer settings from a fitted model or
posterior summaries. A \`NULL\` field represents a setting not recorded
for that result.

## Usage

``` r
make_neurogam_metadata(
  analysis_type,
  time_id,
  space_id = NULL,
  sensor_id = NULL,
  participant_id = NULL,
  outcome_id = NULL,
  predictor_id = NULL,
  multilevel = NULL,
  threshold = NULL,
  threshold_type = NULL,
  null_value = NULL,
  rope = NULL,
  baseline_window = NULL,
  credible_interval = NULL,
  odds_mode = NULL,
  contrast = NULL,
  contrast_levels = NULL,
  adjacency_source = NULL,
  participant_clusters = NULL,
  n_post_samples = NULL
)
```

## Arguments

- analysis_type:

  One of \`"temporal"\`, \`"temporal_generalization"\`, or
  \`"spatiotemporal"\`.

- time_id:

  Names of time coordinates in the resulting predictions.

- space_id:

  Names of projected spatial coordinates, for spatiotemporal results
  only.

- sensor_id:

  Name of the sensor identifier column, for spatiotemporal results only.

- participant_id, outcome_id, predictor_id:

  Original input column names.

- multilevel:

  Model aggregation mode.

- threshold, threshold_type:

  Threshold used for cluster detection.

- null_value:

  Null reference for the effect or contrast.

- rope:

  Fixed ROPE half-width, when one was specified.

- baseline_window:

  Baseline interval defining a data-dependent ROPE, when applicable. The
  resulting numerical ROPE width is not stored here.

- credible_interval:

  Credible interval width.

- odds_mode:

  Posterior-odds strategy actually used by the analysis.

- contrast:

  Description of the estimated effect or difference.

- contrast_levels:

  Explicit condition labels, when available.

- adjacency_source:

  One of \`"auto_knn"\`, \`"provided"\`, or \`NULL\`.

- participant_clusters:

  Whether participant-level clusters were computed.

- n_post_samples:

  Number of posterior draws requested or used.

## Value

A named metadata list.
