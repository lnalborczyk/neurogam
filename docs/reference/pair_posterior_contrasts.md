# Create paired A/B posterior contrasts in long format

Create paired A/B posterior contrasts in long format

## Usage

``` r
pair_posterior_contrasts(
  post_draws,
  coordinates,
  predictor_type = c("categorical", "continuous")
)
```

## Arguments

- post_draws:

  Long predictions with \`.draw\`, \`.epred\`, and \`predictor\`.

- coordinates:

  Time or spatial coordinates; participant optional.

- predictor_type:

  \`categorical\` or \`continuous\`.

## Value

A data frame with \`epred_diff = level2 - level1\` by draw/location.
Categorical factor levels determine order; character values are sorted;
numeric continuous levels are sorted numerically after rounding.
