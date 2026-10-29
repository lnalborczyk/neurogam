# Validate a neurogam analysis result

Validate the common result contract without requiring a fitted Stan
model. This checks structure, not the statistical correctness of the
estimates. Existing concrete S3 classes and legacy fields are not
changed.

## Usage

``` r
validate_neurogam_result(x, check_membership = FALSE)
```

## Arguments

- x:

  A result returned by \`testing_through_time()\` or
  \`analyse_bgam_spatiotemporal()\`.

- check_membership:

  If \`TRUE\`, also verify that cluster membership counts agree with
  \`cluster_summary\$n_points\`. Cluster identity is defined by sign and
  ID, and additionally by participant when applicable.

## Value

Invisibly returns \`TRUE\` or raises an informative error.
