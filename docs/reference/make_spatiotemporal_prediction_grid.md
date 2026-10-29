# Build paired EEG sensor-time prediction grids

Build paired EEG sensor-time prediction grids

## Usage

``` r
make_spatiotemporal_prediction_grid(summary_data, contrast_levels = NULL)
```

## Arguments

- summary_data:

  Output of \`make_spatiotemporal_summary_data()\`.

- contrast_levels:

  Optional two-element character vector defining A - B.

## Value

List with \`a\`, \`b\` grids and a readable contrast label.
