# Select posterior draw indices for prediction

Select posterior draw indices for prediction

## Usage

``` r
select_prediction_draws(n_available, n_requested = NULL)
```

## Arguments

- n_available:

  Number of draws in a fitted model.

- n_requested:

  NULL for all draws, otherwise a positive integer.

## Value

NULL when all draws are requested, otherwise integer draw indices.
