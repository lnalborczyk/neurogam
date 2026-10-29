# Plot spatiotemporal neurogam results with a unified interface

S3 method for 3D sensor-by-time results. Delegates to the existing
\[plot_spatiotemporal_topography()\] and \[plot_sensor_trajectories()\]
implementations; no predictions or clusters are recomputed.

## Usage

``` r
# S3 method for class 'clusters_results_st'
plot(x, type = "topography", observed_data = NULL, observed_data_b = NULL, ...)
```

## Arguments

- x:

  A \`clusters_results_st\` object.

- type:

  One of \`"topography"\` (default), \`"trajectories"\`,
  \`"predictions"\`, or \`"clusters"\`. The latter two use trajectories
  with cluster highlighting switched off or on, respectively.

- observed_data:

  Original long-form EEG observations required by the existing
  spatiotemporal plotting helpers.

- observed_data_b:

  Optional observations for the second condition.

- ...:

  Named parameters passed to the selected specialized plotter.

## Value

A \`ggplot\` object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot(results_st, type = "topography", observed_data = eeg_data2)
plot(results_st, type = "trajectories", observed_data = eeg_data2,
     layout = "facets", cluster_style = "segments")
} # }
```
