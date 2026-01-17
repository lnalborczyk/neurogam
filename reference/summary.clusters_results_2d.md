# Summary method for 2D `clusters_results` objects

Produces a detailed textual summary of a time-resolved Bayesian GAMM
analysis, including model metadata, the number of clusters detected, and
descriptive statistics of cluster durations. A rounded cluster table is
also printed.

## Usage

``` r
# S3 method for class 'clusters_results_2d'
summary(object, digits = 3, ...)
```

## Arguments

- object:

  An object of class `"clusters_results"` created by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).

- digits:

  Integer; number of decimal places used when printing numeric values
  (default: `3`).

- ...:

  Additional arguments (currently ignored).

## Value

The object `object`, returned invisibly.
