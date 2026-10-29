# Find connected threshold-exceeding clusters on arbitrary coordinates

Find connected threshold-exceeding clusters on arbitrary coordinates

## Usage

``` r
find_clusters_nd(
  data,
  coordinates = "time",
  value = "value",
  threshold = 10,
  threshold_type = c("both", "above", "below"),
  group = NULL,
  sensor_id = NULL,
  sensor_adjacency = NULL,
  max_gap = NULL
)
```

## Arguments

- data:

  Data frame containing unique coordinate combinations and a statistic.

- coordinates:

  Character vector of regular-grid coordinate names. For sensor-time
  data, normally \`"time"\` (sensor topology is separate).

- value:

  Column containing positive posterior odds (default \`"value"\`).

- threshold:

  Posterior-odds threshold, strictly greater than one.

- threshold_type:

  \`"both"\`, \`"above"\`, or \`"below"\`.

- group:

  Optional grouping column; connected components never cross groups.

- sensor_id:

  Optional column with sensor identifiers.

- sensor_adjacency:

  Named logical/numeric matrix with sensor row/column names, or
  two-column data frame with sensor pairs.

- max_gap:

  Optional named list of maximum allowed differences for numeric grid
  coordinates (e.g. \`list(time = 0.004)\`). By default adjacent
  observed coordinate levels connect, even when physically far apart.

## Value

List with \`clusters\` (one row per component), \`membership\` (one row
per threshold-exceeding sample), and \`edges\` (full-domain graph
edges).

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.
