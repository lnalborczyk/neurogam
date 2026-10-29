# Construct legacy-compatible temporal results

Internal constructor. Maintains the existing field order and S3 class
names.

## Usage

``` r
assemble_temporal_results(
  clustering,
  predictions,
  model,
  summary_data,
  multilevel,
  is_2d_time,
  metadata = NULL
)
```

## Arguments

- clustering:

  Output from \`assemble_temporal_clusters()\`.

- predictions:

  Posterior summary data frame.

- model:

  Fitted brms model or test fixture.

- summary_data:

  Data used for fitting.

- multilevel:

  Multilevel mode.

- is_2d_time:

  Whether this is a 2D temporal analysis.

## Value

A \`clusters_results_1d\` or \`clusters_results_2d\` object.
