# Prepare a common report for neurogam analysis results

Internal reporting helper used by the 1D, 2D, and spatiotemporal print
and summary methods. Operates only on the standardized result tables and
P2.2 metadata; never refits a model or recomputes posterior quantities.

## Usage

``` r
neurogam_report_data(x)
```

## Arguments

- x:

  A recognized neurogam result.

## Value

A list containing \`metadata\`, \`clusters\`, \`n_clusters\`,
\`model_class\`, and \`n_draws\`.
