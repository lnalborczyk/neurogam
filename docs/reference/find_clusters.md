# Find contiguous clusters in a time series

Identify contiguous clusters of time points in a time series where a
variable exceeds (or falls below) a given threshold. Clusters are
defined as consecutive time points where the thresholding condition
holds.

## Usage

``` r
find_clusters(data, threshold = 10, above_threshold = TRUE)
```

## Arguments

- data:

  A data frame containing at least two columns named `"time"` and
  `"value"`. Additional columns are ignored.

- threshold:

  Numeric scalar; threshold used to define clusters. When
  `above_threshold = TRUE`, clusters are defined where
  `value >= threshold`. When `above_threshold = FALSE`, clusters are
  defined where `value <= threshold`.

- above_threshold:

  Logical scalar; if `TRUE` (default), clusters are formed where
  `value >= threshold`. If `FALSE`, clusters are formed where
  `value <= threshold`.

## Value

A data frame with one row per detected cluster and columns:

- `cluster_id`: integer cluster index (starting at 1);

- `cluster_onset`: time of the first point in the cluster;

- `cluster_offset`: time of the last point in the cluster;

- `n_points`: number of time points in the cluster.

If no clusters are found, an empty data frame with these columns is
returned.

## Details

The function assumes that the `time` variable is numeric and that
consecutive rows correspond to consecutive time points. Internally, the
data are first arranged by `time` and any rows with `NA` in `time` or
`value` are removed.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(666)
df <- data.frame(
  time  = seq(0, 1, length.out = 100),
  value = c(rnorm(40, 0, 1), rnorm(20, 5, 1), rnorm(40, 0, 1))
  )

# find clusters where value >= 3
cl <- find_clusters(df, threshold = 3, above_threshold = TRUE)
cl
} # }
```
