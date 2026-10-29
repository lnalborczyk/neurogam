# Internal dispatch classification for neurogam analysis dimensions

This function routes analyses without altering any data or model
settings. Validation of 1D/2D and spatiotemporal data stays in the
respective analysis implementations.

## Usage

``` r
classify_bgam_analysis(space_id = NULL)
```

## Arguments

- space_id:

  Optional spatial-coordinate column names.

## Value

Either \`"temporal"\` or \`"spatiotemporal"\`.
