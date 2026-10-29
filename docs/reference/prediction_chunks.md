# Split observation positions into consecutive prediction chunks

Split observation positions into consecutive prediction chunks

## Usage

``` r
prediction_chunks(n, chunk_size = 256L)
```

## Arguments

- n:

  Number of prediction locations.

- chunk_size:

  Maximum number of locations in each chunk.

## Value

List of integer position vectors.
