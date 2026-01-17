# Plot time-resolved GAM results (1D clusters)

Visualises the fitted posterior effect through time for a
`clusters_results_1d` object. The plot shows:

- the posterior mean (`x$predictions$post_prob`) as a line,

- its uncertainty interval (`x$predictions$lower`–`x$predictions$upper`)
  as a ribbon,

- and detected significant clusters (`x$clusters`) as thick horizontal
  segments.

## Usage

``` r
# S3 method for class 'clusters_results_1d'
plot(
  x,
  null_value = 0,
  clusters_y = -Inf,
  clusters_colour = "black",
  lineend = "butt",
  theme = ggplot2::theme_bw(),
  ...
)
```

## Arguments

- x:

  A `clusters_results_1d` object (typically returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md)
  for 1D time).

- null_value:

  Numeric. Reference value shown as a horizontal dashed line. Defaults
  to `0`.

- clusters_y:

  Numeric. Vertical position used to draw the cluster segments. Defaults
  to `-Inf` (draws at the bottom of the panel).

- clusters_colour:

  Character. Colour used for the prediction line/ribbon and the cluster
  segments. Defaults to `"black"`.

- lineend:

  Character. Line end style for cluster segments, passed to
  [`geom_segment`](https://ggplot2.tidyverse.org/reference/geom_segment.html)
  (e.g., `"butt"`, `"round"`).

- theme:

  A ggplot2 theme object. Defaults to
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- ...:

  Currently unused. Included for S3 compatibility.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

When possible (depending on `x$multilevel` and the structure of
`x$model$data`), the method also reconstructs and overlays an empirical
time course (e.g., averaged outcome or difference between two predictor
levels).

If clusters are available at the participant level (i.e., a
`participant` column exists in `x$clusters`), the plot is facetted by
participant.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- testing_through_time(...)
plot(res)
} # }
```
