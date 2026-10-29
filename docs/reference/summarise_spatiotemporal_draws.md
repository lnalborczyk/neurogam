# Summarise posterior draws at arbitrary locations

Summarises a matrix of paired posterior contrasts or one-sample
expectations (draws in rows, observation locations in columns).
Posterior odds are directional and can use a symmetric ROPE.

## Usage

``` r
summarise_spatiotemporal_draws(
  draws,
  null_value = 0,
  rope = 0,
  credible_interval = 0.95,
  odds_mode = "directional_rope"
)
```

## Arguments

- draws:

  Numeric matrix, posterior draws by locations.

- null_value:

  Null value against which draws are assessed.

- rope:

  Non-negative half-width of a ROPE around \`null_value\`.

- credible_interval:

  Width of equal-tailed credible intervals.

- odds_mode:

  Posterior-odds inference strategy. The default,
  \`"directional_rope"\`, computes independent continuous posterior odds
  for effects above and below the ROPE. These odds are thresholded
  separately to identify positive and negative clusters.
  \`"legacy_one_sided"\` reproduces the historical capped one-sided
  odds. \`"legacy_directional"\` and \`"spatiotemporal"\` are retained
  for backward compatibility.

## Value

Data frame with posterior estimates, probabilities, and odds.
