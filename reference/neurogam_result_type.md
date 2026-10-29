# Identify a neurogam analysis result

Internal type detection shared by the read-only result accessors. The
existing concrete S3 classes remain authoritative for dispatch.

## Usage

``` r
neurogam_result_type(x)
```

## Arguments

- x:

  A neurogam result object.

## Value

One of \`"temporal"\`, \`"temporal_generalization"\`, or
\`"spatiotemporal"\`.
