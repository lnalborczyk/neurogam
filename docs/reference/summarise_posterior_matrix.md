# Summarise posterior draws at arbitrary prediction locations

Each column is one prediction location, and each row is a paired
posterior draw. This low-level implementation does not know about EEG
coordinates, sensors, conditions, or time. Those are handled by
adapters.

## Usage

``` r
summarise_posterior_matrix(
  draws,
  null_value = 0,
  rope = 0,
  credible_interval = 0.95,
  odds_mode = c("directional_rope", "legacy_one_sided", "legacy_directional",
    "spatiotemporal"),
  odds_cap = NULL
)
```

## Arguments

- draws:

  Finite numeric matrix, posterior draws by locations.

- null_value:

  Point-null value. For \`legacy_directional\`, the null is fixed at
  zero and \`rope\` is the half-width around it.

- rope:

  Non-negative half-width of the symmetric ROPE.

- credible_interval:

  Width of the equal-tailed credible interval.

- odds_mode:

  One of \`legacy_one_sided\`, \`legacy_directional\`, or
  \`spatiotemporal\`.

- odds_cap:

  Finite-sample cap for legacy odds. Defaults to the number of posterior
  draws supplied; ignored for \`spatiotemporal\`.

## Value

Data frame with one row per prediction location.

## Details

The \`legacy_one_sided\` convention is the default across neurogam.
Other conventions remain available explicitly for reproducibility.
Changing the default is an intentional inference-policy change.
