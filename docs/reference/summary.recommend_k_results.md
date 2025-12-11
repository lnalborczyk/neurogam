# Summary method for recommend_k_results objects

Provides a concise summary of k-tuning results, including the range of
tested `k` values, the effective complexity measures (e.g., `p_loo` or
`p_waic`), and the recommended `k` obtained via knee detection.

## Usage

``` r
# S3 method for class 'recommend_k_results'
summary(object, digits = 3, ...)
```

## Arguments

- object:

  An object of class `"recommend_k_results"`.

- digits:

  Integer; number of digits to display for numeric summaries.

- ...:

  Further arguments passed to or from other methods (currently ignored).

## Value

The input object `object`, invisibly.
