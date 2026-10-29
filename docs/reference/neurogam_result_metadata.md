# Access standardized neurogam result metadata

Derives consistent metadata for temporal (1D), cross-temporal (2D), and
sensor-by-time results without mutating the original result. Metadata
absent from older results are kept \`NULL\` rather than guessed. In
particular, the posterior-odds convention is not inferred from the model
or result class.

## Usage

``` r
neurogam_result_metadata(x)
```

## Arguments

- x:

  A result returned by \`testing_through_time()\` or
  \`analyse_bgam_spatiotemporal()\`.

## Value

A named list with analysis type, coordinate names, and any recorded
inference or adjacency settings.
