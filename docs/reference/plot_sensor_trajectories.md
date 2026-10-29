# Compare observed and predicted EEG time courses across sensors

Plots observed sensor-wise means (black) against posterior median model
predictions (blue), optionally with credible ribbons. Trajectories can
be placed in approximate scalp positions, sharing one amplitude scale,
or displayed as ordinary facets. Detected clusters can appear as
translucent rectangles or as thicker, opaque portions of the
\*predicted\* trajectories. Cluster segments respect sign, cluster ID
and gaps in the time grid.

## Usage

``` r
plot_sensor_trajectories(
  results,
  observed_data,
  observed_data_b = NULL,
  outcome_id = "voltage",
  sensor_id = "channel",
  time_id = "time",
  spatial_id = c("xproj", "yproj"),
  participant_id = NULL,
  layout = c("scalp", "facets"),
  cluster_style = c("rectangles", "segments", "none"),
  cluster_alpha = 0.18,
  highlighted_linewidth = 1.1,
  max_gap = NULL,
  show_ci = TRUE,
  ncol = 8L,
  free_y = FALSE,
  sensors = NULL,
  sensor_order = NULL,
  mini_width = NULL,
  mini_height = NULL,
  amplitude_limits = NULL,
  show_labels = TRUE,
  head_outline = TRUE
)
```

## Arguments

- results:

  Output from \[analyse_bgam_spatiotemporal()\].

- observed_data:

  Original long-form observations, or condition A.

- observed_data_b:

  Observed condition B for A - B contrasts.

- outcome_id:

  Observed outcome column.

- sensor_id:

  Sensor identifier column.

- time_id:

  Numeric time column.

- spatial_id:

  Two projected sensor coordinates.

- participant_id:

  Optional participant ID for an equal-participant- weighted observed
  average. With \`NULL\`, observations are pooled.

- layout:

  \`"scalp"\` for mini-trajectories at electrode positions or
  \`"facets"\` for one independent panel per sensor.

- cluster_style:

  \`"rectangles"\` for transparent backgrounds, \`"segments"\` for
  highlighted prediction lines, or \`"none"\`.

- cluster_alpha:

  Opacity of cluster rectangle backgrounds.

- highlighted_linewidth:

  Width of cluster-highlighted segments.

- max_gap:

  Maximum time gap for joining cluster members. Defaults to 1.5 times
  the median sampling interval.

- show_ci:

  Plot posterior credible interval ribbons.

- ncol:

  Number of facet columns when \`layout = "facets"\`.

- free_y:

  Use independent amplitude limits in faceted plots.

- sensors:

  Optional vector of sensors to display.

- sensor_order:

  Optional vector specifying the order of facet panels.

- mini_width:

  Width of each scalp mini-trace in projected units; defaults to 80
  percent of median nearest-sensor spacing.

- mini_height:

  Height of each scalp mini-trace in projected units; defaults to 75
  percent of median nearest-sensor spacing.

- amplitude_limits:

  Optional common amplitude range for scalp traces.

- show_labels:

  Display sensor labels in scalp layout.

- head_outline:

  Add a schematic head circle and nose in scalp layout.

## Value

A ggplot object.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_sensor_trajectories(
    results_st, eeg_data2, outcome_id = "voltage",
    sensor_id = "channel", layout = "scalp",
    cluster_style = "segments"
)
} # }
```
