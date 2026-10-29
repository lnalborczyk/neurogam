# Explicit plotting options for temporal results

Called only when the user supplies \`type\` to \`plot()\` for a 1D
result. The historical no-\`type\` plot remains in R/plot.R, unchanged.

## Usage

``` r
neurogam_plot_typed_1d(
  x,
  type,
  null_value = 0,
  clusters_colour = "black",
  theme = theme_neurogam(),
  ...
)
```
