# Plot spatio-temporal EEG data as topographies or 2D surfaces

Create facetted EEG maps across multiple time points from either i) a
brms spatio-temporal model (predicted values) or ii) raw/summarised EEG
data in long format (observed values).

## Usage

``` r
plot_eeg(
  x,
  type = c("topo", "surface"),
  times = NULL,
  n_times = NULL,
  sensors = NULL,
  value_col = "voltage",
  time_col = "time",
  x_col = "xproj",
  y_col = "yproj",
  grid_res = 100,
  head_expand = 1.1,
  re_formula = NULL,
  probs = c(0.025, 0.975),
  ndraws = NULL,
  show_sensors = TRUE,
  sensor_size = 1.5,
  sensor_labels = FALSE,
  sensor_label_col = "channel",
  sensor_label_size = 3,
  sensor_label_repel = TRUE,
  contours = TRUE,
  contour_bins = 10,
  palette = "RdBu",
  facet_nrow = NULL,
  facet_ncol = NULL,
  facet_scales = "fixed",
  facet_label_prefix = "Time: ",
  facet_unit = "s",
  fill_limits = c("global_quantile", "global", "none"),
  limit_quantiles = c(0.01, 0.99)
)
```

## Arguments

- x:

  Either a data frame containing EEG data (raw or summarised) or a
  fitted `brmsfit` object. If a data frame is provided, it must contain
  columns specified by `time_col`, `x_col`, `y_col`, and `value_col`. If
  a `brmsfit` is provided, predictions are computed from `x$data` and
  `x` must include a compatible spatio-temporal smooth (e.g.,
  `gp(time, xproj, yproj, ...)` or `t2(time, xproj, yproj, ...)`).

- type:

  Character; plot type. `"topo"` draws a scalp outline and masks values
  outside the head. `"surface"` draws only the interpolated/predicted
  field.

- times:

  Numeric vector of time points to plot. If `NULL`, all unique times
  found in the input are used (or a subset if `n_times` is provided).

- n_times:

  Optional integer; if not `NULL`, selects `n_times` approximately
  equally spaced time points from `times`.

- sensors:

  Optional data frame of sensor positions used to draw the head outline,
  sensor points, and sensor labels. Must contain columns `xproj` and
  `yproj`. If `NULL`, sensor coordinates are inferred from the unique
  `x_col`/`y_col` pairs in `x` (labels will then be unavailable).

- value_col:

  Character; name of the value column in raw data (e.g., `"voltage"`).

- time_col, x_col, y_col:

  Character; column names in `x` giving time and 2D sensor coordinates
  (projected).

- grid_res:

  Integer; grid resolution used to evaluate predictions and/or
  interpolate data to a regular grid (approximately `grid_res^2` pixels
  per facet).

- head_expand:

  Numeric; expansion factor used to compute the head radius from sensor
  coordinates (values \> 1 add padding).

- re_formula:

  Passed to
  [`brms::fitted`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)
  when `x` is a `brmsfit`. Use `NULL` for default behaviour, or `NA` to
  exclude group-level terms.

- probs:

  Numeric vector of length 2 giving the quantiles returned by
  `brms::fitted(..., probs = probs)`.

- ndraws:

  Optional integer; number of posterior draws used by `brms::fitted`.

- show_sensors:

  Logical; if `TRUE` and `type = "topo"`, plot sensor points.

- sensor_size:

  Numeric; point size for sensors.

- sensor_labels:

  Logical; if `TRUE`, add sensor labels (requires `sensors` to contain
  `sensor_label_col`).

- sensor_label_col:

  Character; column name in `sensors` providing sensor names.

- sensor_label_size:

  Numeric; text size for sensor labels.

- sensor_label_repel:

  Logical; if `TRUE`, use ggrepel for labels.

- contours:

  Logical; if `TRUE`, add contour lines.

- contour_bins:

  Integer; number of contour bins.

- palette:

  Character; palette name passed to
  [`scale_fill_distiller`](https://ggplot2.tidyverse.org/reference/scale_brewer.html).

- facet_nrow, facet_ncol:

  Integers; layout for
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

- facet_scales:

  Character; facet scaling, passed to `facet_wrap(scales = ...)`.

- facet_label_prefix:

  Character; prefix used when generating facet labels (e.g.,
  `"Time: "`).

- facet_unit:

  Character; unit suffix used when generating facet labels (e.g.,
  `"s"`).

- fill_limits:

  Character; how to set colour scale limits. `"global_quantile"` uses
  `limit_quantiles` to define robust symmetric limits, `"global"` uses
  the full range, and `"none"` leaves limits to ggplot.

- limit_quantiles:

  Numeric vector of length 2; quantiles used when
  `fill_limits = "global_quantile"`.

## Value

A `ggplot` object.

## Details

The function can render either a scalp topography (`type = "topo"`) with
a head outline and optional sensor markers/labels, or a simple 2D
surface plot (`type = "surface"`) without a head outline.

For `type = "topo"`, values are masked outside the head radius computed
from sensor coordinates and a head outline (circle + nose + ears) is
drawn.

For raw data, values are interpolated to a regular grid per time point
using `st_interp_to_grid()`, which relies on akima and can optionally
fill missing grid cells using nearest-neighbour (FNN).

For `brmsfit` inputs, predictions are obtained via `brms::fitted()` on a
regular grid. When `type = "topo"`, missing values inside the head can
be filled using `st_fill_head()` (nearest-neighbour).

This function relies on internal helpers such as
[`st_take_n_times()`](https://lnalborczyk.github.io/neurogam/reference/st_take_n_times.md),
`st_make_grid()`, `st_head_radius()`, `st_head_outline()`,
`st_build_time_labels()`, `st_order_time_factor()`,
`st_compute_limits()`, `st_predict_brms()`, `st_fill_head()`, and
`st_interp_to_grid()`.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>

## Examples

``` r
if (FALSE) { # \dontrun{
# Summarised EEG data (long format)
plot_eeg(
  eeg_data_summary,
  type = "topo",
  sensors = sensors,
  times = c(0, 0.1, 0.2, 0.3),
  grid_res = 80,
  contours = FALSE,
  facet_nrow = 2
  )

# brms model predictions
plot_eeg(
  spatio_temporal_gam,
  type = "topo",
  sensors = sensors,
  times = c(0, 0.1, 0.2, 0.3),
  ndraws = 200
  )
} # }
```
