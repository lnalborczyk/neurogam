# Find contiguous clusters in a time series

Identify contiguous clusters of time points in a time series where a
variable exceeds a positive threshold or falls below a negative
threshold. Clusters are defined as consecutive time points satisfying
the thresholding condition.

## Usage

``` r
find_clusters(
  data,
  threshold = 10,
  group = NULL,
  above_threshold = TRUE,
  both_signs = TRUE
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
  and negative clusters where `value <= -threshold`. Must be
  non-negative when `both_signs = TRUE`.

- group:

  Optional grouping column name (character scalar) used to find clusters
  independently within each group level (e.g., `"participant"`). Set to
  `NULL` (default) to ignore grouping.

- above_threshold:

  Logical scalar used only when `both_signs = FALSE`. If `TRUE`
  (default), clusters are formed where `value >= threshold`; if `FALSE`,
  clusters are formed where `value <= threshold`.

- both_signs:

  Logical scalar indicating whether to detect both positive and negative
  clusters in a single call (default: `TRUE`). When `TRUE`,
  `above_threshold` is ignored.

## Value

A data frame with one row per detected cluster and columns:

- `cluster_id`: integer cluster index (starting at 1). When `group` is
  provided, `cluster_id` restarts at 1 within each group level;

- `cluster_onset`: time of the first point in the cluster;

- `cluster_offset`: time of the last point in the cluster;

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

## Examples

``` r
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
find_clusters(data = df, threshold = 3)
#>   cluster_id cluster_onset cluster_offset n_points     sign
#> 1          1     0.3030303      0.3636364        7 positive
#> 2          2     0.3838384      0.4444444        7 positive
#> 3          3     0.4646465      0.4949495        4 positive
#> 4          1     0.5050505      0.5050505        1 negative
#> 5          2     0.5252525      0.5858586        7 negative
#> 6          3     0.6060606      0.6262626        3 negative
#> 7          4     0.6464646      0.6767677        4 negative
#> 8          5     0.6969697      0.6969697        1 negative

# One-sided detection (positive only)
find_clusters(data = df, threshold = 3, both_signs = FALSE, above_threshold = TRUE)
#>   cluster_id cluster_onset cluster_offset n_points     sign
#> 1          1     0.3030303      0.3636364        7 positive
#> 2          2     0.3838384      0.4444444        7 positive
#> 3          3     0.4646465      0.4949495        4 positive

# Grouped example (e.g., per participant)
df_g <- rbind(
  transform(df, participant = "P01"),
  transform(df, participant = "P02")
  )

find_clusters(
  data = df_g,
  threshold = 3,
  group = "participant"
  )
#>    participant cluster_id cluster_onset cluster_offset n_points     sign
#> 1          P01          1     0.3030303      0.3636364        7 positive
#> 2          P01          2     0.3838384      0.4444444        7 positive
#> 3          P01          3     0.4646465      0.4949495        4 positive
#> 4          P02          1     0.3030303      0.3636364        7 positive
#> 5          P02          2     0.3838384      0.4444444        7 positive
#> 6          P02          3     0.4646465      0.4949495        4 positive
#> 7          P01          1     0.5050505      0.5050505        1 negative
#> 8          P01          2     0.5252525      0.5858586        7 negative
#> 9          P01          3     0.6060606      0.6262626        3 negative
#> 10         P01          4     0.6464646      0.6767677        4 negative
#> 11         P01          5     0.6969697      0.6969697        1 negative
#> 12         P02          1     0.5050505      0.5050505        1 negative
#> 13         P02          2     0.5252525      0.5858586        7 negative
#> 14         P02          3     0.6060606      0.6262626        3 negative
#> 15         P02          4     0.6464646      0.6767677        4 negative
#> 16         P02          5     0.6969697      0.6969697        1 negative
```
