# Plot time-generalisation GAM results (2D clusters)

Visualises 2D (time-by-time) posterior predictions from a
`clusters_results_2d` object as a two-panel figure:

1.  summary empirical data

2.  summary posterior predictions

## Usage

``` r
# S3 method for class 'clusters_results_2d'
plot(
  x,
  theme = ggplot2::theme_bw(),
  palette = "vik",
  midpoint = 0.5,
  axes_labels = c("Training time (s)", "Testing time (s)"),
  cluster_outline_colour = "black",
  ...
)
```

## Arguments

- x:

  A `clusters_results_2d` object (typically returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md)
  for 2D time / time-generalisation).

- theme:

  A ggplot2 theme object added to each panel. Defaults to
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- palette:

  Character, a scico palette, see scico::scico_palette_show().

- midpoint:

  Numeric, the midpoint for diverging color palettes.

- axes_labels:

  Character vector, the x- and y-axis labels.

- cluster_outline_colour:

  Character. Colour used to outline cluster cells. Defaults to
  `"white"`.

- ...:

  Currently unused. Included for S3 compatibility.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
(a patchwork composition of two ggplots).

## Details

Significant clusters stored in `x$clusters` are assumed to be
**pointwise** (one row per time-by-time cell). Cluster locations are
overlaid as contour (which includes interpolation).

This method relies on patchwork to combine the two panels using `+`.

## Examples

``` r
if (FALSE) { # \dontrun{
res2d <- testing_through_time(...)
plot(res2d)
} # }
```
