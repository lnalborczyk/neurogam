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
  threshold_type = c("both", "above", "below")
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
find_clusters(data = df, threshold = 3, threshold_type = "both")
#>    id      onset     offset n_points     sign
#> 1   1 0.30303030 0.36363636        7 positive
#> 2   2 0.38383838 0.44444444        7 positive
#> 3   3 0.46464646 0.49494949        4 positive
#> 4   1 0.02020202 0.02020202        1 negative
#> 5   2 0.04040404 0.04040404        1 negative
#> 6   3 0.06060606 0.09090909        4 negative
#> 7   4 0.11111111 0.11111111        1 negative
#> 8   5 0.13131313 0.15151515        3 negative
#> 9   6 0.18181818 0.18181818        1 negative
#> 10  7 0.20202020 0.21212121        2 negative
#> 11  8 0.23232323 0.29292929        7 negative
#> 12  9 0.50505051 0.71717172       22 negative
#> 13 10 0.73737374 0.79797980        7 negative
#> 14 11 0.82828283 0.84848485        3 negative
#> 15 12 0.88888889 0.94949495        7 negative
#> 16 13 0.97979798 1.00000000        3 negative

# One-sided detection (positive only)
find_clusters(data = df, threshold = 3, threshold_type = "above")
#>   id     onset    offset n_points     sign
#> 1  1 0.3030303 0.3636364        7 positive
#> 2  2 0.3838384 0.4444444        7 positive
#> 3  3 0.4646465 0.4949495        4 positive

# One-sided detection (negative only)
find_clusters(data = df, threshold = 3, threshold_type = "below")
#>    id      onset     offset n_points     sign
#> 1   1 0.02020202 0.02020202        1 negative
#> 2   2 0.04040404 0.04040404        1 negative
#> 3   3 0.06060606 0.09090909        4 negative
#> 4   4 0.11111111 0.11111111        1 negative
#> 5   5 0.13131313 0.15151515        3 negative
#> 6   6 0.18181818 0.18181818        1 negative
#> 7   7 0.20202020 0.21212121        2 negative
#> 8   8 0.23232323 0.29292929        7 negative
#> 9   9 0.50505051 0.71717172       22 negative
#> 10 10 0.73737374 0.79797980        7 negative
#> 11 11 0.82828283 0.84848485        3 negative
#> 12 12 0.88888889 0.94949495        7 negative
#> 13 13 0.97979798 1.00000000        3 negative

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
#>    participant id      onset     offset n_points     sign
#> 1          P01  1 0.30303030 0.36363636        7 positive
#> 2          P01  2 0.38383838 0.44444444        7 positive
#> 3          P01  3 0.46464646 0.49494949        4 positive
#> 4          P01  1 0.02020202 0.02020202        1 negative
#> 5          P01  2 0.04040404 0.04040404        1 negative
#> 6          P01  3 0.06060606 0.09090909        4 negative
#> 7          P01  4 0.11111111 0.11111111        1 negative
#> 8          P01  5 0.13131313 0.15151515        3 negative
#> 9          P01  6 0.18181818 0.18181818        1 negative
#> 10         P01  7 0.20202020 0.21212121        2 negative
#> 11         P01  8 0.23232323 0.29292929        7 negative
#> 12         P01  9 0.50505051 0.71717172       22 negative
#> 13         P01 10 0.73737374 0.79797980        7 negative
#> 14         P01 11 0.82828283 0.84848485        3 negative
#> 15         P01 12 0.88888889 0.94949495        7 negative
#> 16         P01 13 0.97979798 1.00000000        3 negative
#> 17         P02  1 0.30303030 0.36363636        7 positive
#> 18         P02  2 0.38383838 0.44444444        7 positive
#> 19         P02  3 0.46464646 0.49494949        4 positive
#> 20         P02  1 0.02020202 0.02020202        1 negative
#> 21         P02  2 0.04040404 0.04040404        1 negative
#> 22         P02  3 0.06060606 0.09090909        4 negative
#> 23         P02  4 0.11111111 0.11111111        1 negative
#> 24         P02  5 0.13131313 0.15151515        3 negative
#> 25         P02  6 0.18181818 0.18181818        1 negative
#> 26         P02  7 0.20202020 0.21212121        2 negative
#> 27         P02  8 0.23232323 0.29292929        7 negative
#> 28         P02  9 0.50505051 0.71717172       22 negative
#> 29         P02 10 0.73737374 0.79797980        7 negative
#> 30         P02 11 0.82828283 0.84848485        3 negative
#> 31         P02 12 0.88888889 0.94949495        7 negative
#> 32         P02 13 0.97979798 1.00000000        3 negative
```
