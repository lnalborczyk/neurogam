# Summary method for `clusters_results` objects

Produces a detailed textual summary of a time-resolved Bayesian GAMM
analysis, including model metadata, the number of clusters detected, and
descriptive statistics of cluster durations. A rounded cluster table is
also printed.

## Usage

``` r
# S3 method for class 'clusters_results'
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

## Details

The summary prints:

- the model type used (`"full"`, `"summary"`, or `"group"`);

- the class of the underlying brms model and the number of posterior
  draws;

- the number of clusters detected by the posterior odds threshold;

- descriptive statistics for cluster durations (minimum, maximum, mean,
  median, and total duration);

- a neatly formatted table listing each cluster's onset, offset, and
  duration.

If no clusters were detected, the function prints a message and returns
invisibly.

## See also

[`print.clusters_results`](https://lnalborczyk.github.io/neurogam/reference/print.clusters_results.md),
[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md)
