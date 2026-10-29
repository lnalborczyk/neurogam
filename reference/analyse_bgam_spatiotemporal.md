# Identify posterior-odds clusters in sensor-by-time BGAMM predictions

Uses an existing univariate \`brmsfit\` object and explicit prediction
grids. When sensor adjacency is omitted, an approximate adjacency graph
is inferred from projected sensor coordinates. An explicitly supplied
physical sensor topology is preferable for confirmatory analyses.

## Usage

``` r
analyse_bgam_spatiotemporal(
  model,
  newdata,
  sensor_adjacency = NULL,
  newdata_b = NULL,
  sensor_id = "sensor",
  time_id = "time",
  spatial_id = c("xproj", "yproj"),
  null_value = 0,
  rope = 0,
  threshold = 10,
  threshold_type = c("both", "above", "below"),
  credible_interval = 0.95,
  draw_ids = NULL,
  chunk_size = 256L,
  re_formula = NA,
  max_gap = NULL,
  adjacency_k = 4L,
  adjacency_distance_factor = 1.8,
  odds_mode = "directional_rope"
)
```

## Arguments

- model:

  Fitted univariate \`brmsfit\` model.

- newdata:

  Data frame containing exactly one row per sensor and time. Additional
  model predictor columns, including grouping factors, must be present
  as required by \`brms::posterior_epred()\`.

- sensor_adjacency:

  Named adjacency matrix, two-column sensor edge list, or \`NULL\`
  (default) to infer adjacency from projected sensor positions.

- newdata_b:

  Optional second condition grid with exactly the same rows and
  sensor/time coordinates as \`newdata\`. Inference then concerns A - B.

- sensor_id:

  Name of the electrode identifier column.

- time_id:

  Name of the time coordinate column.

- spatial_id:

  Names of two projected coordinate columns (for plotting).

- null_value:

  Null value for the expected outcome or contrast.

- rope:

  Nonnegative half-width of a fixed ROPE around \`null_value\`.

- threshold:

  Posterior odds threshold, strictly greater than one.

- threshold_type:

  Whether to detect above/below/both directions.

- credible_interval:

  Width of equal-tailed posterior intervals.

- draw_ids:

  Optional posterior draw indices. When omitted, use all draws.
  Explicitly supplying indices ensures pairing across conditions and
  chunks.

- chunk_size:

  Number of sensor-time locations processed at once.

- re_formula:

  \`NA\` for population-level predictions (default), \`NULL\` to include
  fitted group-level deviations.

- max_gap:

  Optional \`list(time = ...)\` passed to \`find_clusters_nd()\`.

- adjacency_k:

  Number of nearest neighbors considered when adjacency is generated
  automatically.

- adjacency_distance_factor:

  Maximum auto-generated edge length relative to the median
  nearest-neighbor distance. Set to \`NULL\` to disable distance
  pruning.

- odds_mode:

  Posterior-odds inference strategy. The default,
  \`"directional_rope"\`, computes independent continuous posterior odds
  for effects above and below the ROPE. These odds are thresholded
  separately to identify positive and negative clusters.
  \`"legacy_one_sided"\` reproduces the historical capped one-sided
  odds. \`"legacy_directional"\` and \`"spatiotemporal"\` are retained
  for backward compatibility.

## Value

List with \`predictions\`, \`cluster_summary\`, \`cluster_membership\`,
\`model\`, \`sensor_adjacency\`, \`adjacency_source\`, and metadata.

## See also

\[make_sensor_adjacency()\], \[plot_sensor_adjacency()\]
