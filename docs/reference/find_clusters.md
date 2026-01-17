# Find contiguous clusters in a time series

Identify contiguous clusters of time points in a time series where a
variable exceeds a positive threshold, falls below a negative threshold,
or both. Clusters are defined as consecutive time points satisfying the
thresholding condition.

## Usage

``` r
find_clusters(
  data,
  threshold = 10,
  group = NULL,
  threshold_type = c("both", "above", "below"),
  time_id = "time"
)
```

## Arguments

- data:

  A data frame containing at least two columns named `"time"` and
  `"value"`. If `group` is not `NULL`, `data` must also contain that
  grouping column.

- threshold:

  Numeric scalar specifying the (positive) threshold used to define
  clusters. Positive clusters are defined where `value >= threshold`,
  and negative clusters where `value <= 1/threshold`. Must be
  non-negative when `threshold_type = "both"`.

- group:

  Optional grouping column name (character scalar) used to find clusters
  independently within each group level (e.g., `"participant"`). Set to
  `NULL` (default) to ignore grouping.

- threshold_type:

  Character scalar controlling which clusters are detected. Must be one
  of `"above"`, `"below"`, or `"both"` (default). When `"above"`,
  clusters are formed where `value >= threshold`. When `"below"`,
  clusters are formed where `value <= 1/threshold`. When `"both"`, both
  types are detected and the returned data include a `sign` column.

- time_id:

  Character; name of the column(s) in `data` containing time information
  (e.g., in seconds or samples).

## Value

A data frame with one row per detected cluster and columns:

- `id`: integer cluster index (starting at 1). When `group` is provided,
  `id` restarts at 1 within each group level;

- `onset`: time of the first point in the cluster;

- `offset`: time of the last point in the cluster;

- `n_points`: number of time points in the cluster;

- `sign`: character indicating cluster sign (`"positive"` or
  `"negative"`).

If `group` is not `NULL`, the returned data frame also contains the
grouping column (named as in `group`).

If no clusters are found, an empty data frame with the same column
structure is returned.

## Details

By default, the function detects **both positive and negative clusters**
in a single call and returns a column indicating the cluster sign.

If a grouping variable is provided (e.g., `"participant"`), clusters are
detected independently within each group level.

The function assumes that the `time` variable is numeric and that
consecutive rows correspond to consecutive time points (within each
group if grouping is used). Internally, the data are:

1.  filtered to remove rows with missing values;

2.  arranged by `time` (and by `group` then `time`, if used);

3.  thresholded to identify positive and/or negative excursions;

4.  segmented into runs of consecutive threshold-exceeding values, which
    define clusters.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(666)
df <- data.frame(
  time = seq(0, 1, length.out = 100),
  value = c(
    rnorm(30, 0, 1),
    rnorm(20,  4, 1), # positive cluster
    rnorm(20, -4, 1), # negative cluster
    rnorm(30, 0, 1)
    )
  )

# Detect both positive and negative clusters
find_clusters(data = df, threshold = 3, threshold_type = "both")

# One-sided detection (positive only)
find_clusters(data = df, threshold = 3, threshold_type = "above")

# One-sided detection (negative only)
find_clusters(data = df, threshold = 3, threshold_type = "below")

# Grouped example (e.g., per participant)
df_g <- rbind(
  transform(df, participant = "P01"),
  transform(df, participant = "P02")
  )

find_clusters(
  data = df_g,
  threshold = 3,
  group = "participant",
  threshold_type = "both"
  )
} # }
```
