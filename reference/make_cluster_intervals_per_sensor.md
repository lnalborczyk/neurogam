# Construct nonbridging sensor-specific cluster intervals

Splits cluster memberships at changes of sign, cluster ID, absent time
samples, or time gaps exceeding \`max_gap\`.

## Usage

``` r
make_cluster_intervals_per_sensor(data, max_gap = NULL)
```

## Arguments

- data:

  Standardized output of \`prepare_spatiotemporal_plot_data()\`.

- max_gap:

  Largest allowed temporal gap. By default, 1.5 times the median spacing
  in the prediction grid.

## Value

Data frame containing sensor, sign, cluster_id, xmin and xmax.
