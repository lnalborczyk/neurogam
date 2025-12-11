# Print method for `clusters_results` objects

This method provides a concise console representation of the output from
[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
including the number of detected clusters and a compact table
summarising each cluster's onset, offset, and duration. Values are
rounded for readability.

## Usage

``` r
# S3 method for class 'clusters_results'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"clusters_results"` as returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).

- digits:

  Integer; number of decimal places used when printing numeric values
  (default: `3`).

- ...:

  Additional arguments (currently ignored).

## Value

The input object `x`, returned invisibly.

## Details

The printed cluster table includes:

- `cluster_id`: numeric identifier of the cluster;

- `cluster_onset`: estimated temporal onset of the cluster;

- `cluster_offset`: estimated temporal offset of the cluster;

- `duration`: duration of the cluster, computed as
  `cluster_offset - cluster_onset`.

If no clusters exceed the posterior odds threshold, an informative
message is displayed and no table is printed.

## See also

[`summary.clusters_results`](https://lnalborczyk.github.io/neurogam/reference/summary.clusters_results.md),
[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md)
