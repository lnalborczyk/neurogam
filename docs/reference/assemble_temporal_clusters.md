# Detect temporal clusters from posterior summaries

Internal compatibility boundary between posterior inference and
clustering. Preserves the historical 1D/2D cluster output alongside the
standardized cluster summary and membership tables.

## Usage

``` r
assemble_temporal_clusters(
  predictions,
  is_2d_time,
  participant_clusters,
  threshold,
  threshold_type,
  odds_mode = "directional_rope"
)
```

## Arguments

- predictions:

  Data frame containing coordinates and \`prob_ratio\`.

- is_2d_time:

  Whether predictions use \`time1\` and \`time2\` instead of \`time\`.

- participant_clusters:

  Whether predictions are grouped by participant.

- threshold:

  Posterior-odds threshold.

- threshold_type:

  One of \`"both"\`, \`"above"\`, or \`"below"\`.

## Value

A list containing \`clusters\`, \`cluster_summary\`, and
\`cluster_membership\`.
