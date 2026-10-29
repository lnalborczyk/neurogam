# Construct approximate EEG sensor adjacency from projected positions

Builds a symmetric, named logical adjacency matrix using spatial
k-nearest neighbors. To avoid connecting distant sensors, candidate
edges longer than \`distance_factor\` times the median nearest-neighbor
distance are omitted. The result is an \*approximation\* based on 2D
positions; inspect the graph before scientific use.

## Usage

``` r
make_sensor_adjacency(
  data,
  sensor_id = "sensor",
  spatial_id = c("xproj", "yproj"),
  k = 4L,
  distance_factor = 1.8,
  warn_disconnected = TRUE
)
```

## Arguments

- data:

  Data frame containing sensor names and projected locations. Repeated
  rows per sensor (e.g., a sensor-by-time prediction grid) are OK
  provided each sensor has a single consistent position.

- sensor_id:

  Sensor identifier column (default \`"sensor"\`).

- spatial_id:

  Two numeric spatial coordinate columns.

- k:

  Number of nearest neighbors considered per electrode. The resulting
  undirected graph may have more than \`k\` edges per sensor.

- distance_factor:

  Maximum candidate edge length divided by the median nearest-neighbor
  distance. Use \`NULL\` to disable filtering.

- warn_disconnected:

  Issue a warning if the inferred graph has more than one connected
  component.

## Value

Symmetric logical matrix with sensor names on both axes and \`FALSE\` on
the diagonal.
