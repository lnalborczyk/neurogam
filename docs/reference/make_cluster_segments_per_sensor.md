# Construct highlighted prediction-line segments within clusters

Each segment connects only neighboring points sharing the same sensor,
sign and cluster ID; no lines cross missing samples or long time gaps.

## Usage

``` r
make_cluster_segments_per_sensor(data, max_gap = NULL)
```

## Arguments

- data:

  Standardized posterior/observation table.

- max_gap:

  Largest allowed temporal gap.

## Value

Data frame with segment endpoints and sign.
