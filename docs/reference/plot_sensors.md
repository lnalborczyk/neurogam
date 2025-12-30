# Plot EEG sensor positions on a 2D scalp

Visualise the 2D layout of EEG sensors using projected coordinates
(e.g., `xproj` and `yproj`). Optionally draws a head outline, sensor
markers, and/or sensor labels. A subset of sensors can be highlighted
(by name or by row index), with non-highlighted sensors optionally
dimmed.

## Usage

``` r
plot_sensors(
  sensors,
  x_col = "xproj",
  y_col = "yproj",
  label_col = "channel",
  show_head = TRUE,
  head_expand = 1.1,
  show_points = FALSE,
  point_size = 2,
  show_labels = TRUE,
  label_size = 3,
  label_repel = TRUE,
  label_only_highlight = FALSE,
  highlight = NULL,
  highlight_col = "orangered",
  highlight_size = 3,
  highlight_shape = 16,
  highlight_alpha = 1,
  dim_others = TRUE,
  other_alpha = 0.3,
  xlim = NULL,
  ylim = NULL,
  theme = ggplot2::theme_void()
)
```

## Arguments

- sensors:

  A data frame containing sensor coordinates. Must include the columns
  specified by `x_col` and `y_col`. If `show_labels = TRUE` and labels
  are desired, it must also include `label_col`.

- x_col, y_col:

  Character scalars giving the column names in `sensors` corresponding
  to the x and y projected sensor coordinates.

- label_col:

  Character scalar giving the column name in `sensors` used as sensor
  labels (e.g., `"channel"`).

- show_head:

  Logical; if `TRUE`, draw a head outline (circle + nose + ears).

- head_expand:

  Numeric scalar; multiplier applied to the maximum sensor radius when
  computing the head radius. Values \> 1 add padding around sensors.

- show_points:

  Logical; if `TRUE`, draw sensor markers for non-highlighted sensors.

- point_size:

  Numeric; size of non-highlighted sensor markers.

- show_labels:

  Logical; if `TRUE`, add sensor labels.

- label_size:

  Numeric; text size for sensor labels.

- label_repel:

  Logical; if `TRUE`, use ggrepel to reduce label overlap. If `FALSE`,
  uses
  [`geom_text`](https://ggplot2.tidyverse.org/reference/geom_text.html).

- label_only_highlight:

  Logical; if `TRUE` and `highlight` is not `NULL`, show labels only for
  highlighted sensors.

- highlight:

  Optional vector specifying sensors to highlight. If numeric, treated
  as row indices of `sensors`. Otherwise treated as sensor names that
  are matched against `sensors[[label_col]]`.

- highlight_col:

  Character; colour used for highlighted sensors.

- highlight_size:

  Numeric; marker size for highlighted sensors.

- highlight_shape:

  Numeric; point shape for highlighted sensors (as in `ggplot2`).

- highlight_alpha:

  Numeric in \[0, 1\]; alpha transparency for highlighted sensors.

- dim_others:

  Logical; if `TRUE`, apply `other_alpha` to non-highlighted sensor
  markers when `show_points = TRUE`.

- other_alpha:

  Numeric in \[0, 1\]; alpha transparency for non-highlighted sensors
  when `dim_others = TRUE`.

- xlim, ylim:

  Optional numeric vectors of length 2 giving x and y limits. If `NULL`,
  limits are computed automatically from the head radius (if
  `show_head = TRUE`) or from the range of sensor coordinates otherwise.

- theme:

  A [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) object
  modifying the appearance of the plots.

## Value

A `ggplot` object.

## Details

The head outline is computed from the maximum sensor radius (multiplied
by `head_expand`) using `st_head_radius()` and `st_head_outline()`.
These helper functions must be available in the calling environment.

Highlighting works in two modes:

- If `highlight` is numeric, indices are used to mark rows of `sensors`.

- If `highlight` is character, values are matched against
  `sensors[[label_col]]`.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic sensor layout with head outline and labels
plot_sensors(sensors, show_points = TRUE, show_labels = TRUE)

# Highlight a few channels by name
plot_sensors(
  sensors,
  show_points = TRUE,
  show_labels = TRUE,
  highlight = c("CZ", "C1", "C2"),
  label_only_highlight = TRUE
  )

# Highlight by row index
plot_sensors(
  sensors,
  show_points = TRUE,
  highlight = c(10, 25, 42),
  highlight_col = "dodgerblue"
  )
} # }
```
