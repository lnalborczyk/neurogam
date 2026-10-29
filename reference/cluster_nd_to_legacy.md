# Convert multidimensional clustering results to legacy output

Internal adapter retaining the public \`clusters\` data-frame contract.
The standardised cluster summaries and memberships remain available
separately from the output of \`find_clusters_nd()\`.

## Usage

``` r
cluster_nd_to_legacy(result, time_id = "time", group = NULL)
```

## Arguments

- result:

  A list returned by \`find_clusters_nd()\`.

- time_id:

  One or two temporal coordinate column names.

- group:

  An optional grouping column name.

## Value

A legacy-compatible data frame: one row per cluster for 1D, or one row
per cluster member for 2D.
