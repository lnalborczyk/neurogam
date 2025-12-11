# Print method for recommend_k_results objects

Pretty printer for objects of class `"recommend_k_results"` returned by
the k-recommendation procedure (e.g., based on the effective number of
parameters `p_loo` or `p_waic`).

## Usage

``` r
# S3 method for class 'recommend_k_results'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `"recommend_k_results"`.

- digits:

  Integer; number of digits to display for numeric summaries.

- ...:

  Further arguments passed to or from other methods (currently ignored).

## Value

The input object `x`, invisibly.
