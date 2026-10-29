# Create a standardized read-only view of neurogam results

Creates a common interface to predictions and detected clusters without
changing the input object's class, list fields, or legacy \`clusters\`
table. Use this as an intermediate compatibility layer before
introducing unified S3 \`print()\`, \`summary()\`, and \`plot()\`
methods.

## Usage

``` r
as_neurogam_result_view(x, check_membership = FALSE)
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

A list of class \`neurogam_result_view\` containing \`predictions\`,
\`cluster_summary\`, \`cluster_membership\`, \`model\`,
\`summary_data\`, and \`metadata\`. It holds the existing objects, not
recomputed values.
