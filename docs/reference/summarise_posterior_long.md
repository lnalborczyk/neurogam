# Group long-format posterior draws and use the shared inference core

Group long-format posterior draws and use the shared inference core

## Usage

``` r
summarise_posterior_long(data, groups, value = ".epred", ...)
```

## Arguments

- data:

  Long data frame with coordinate columns and one draw value.

- groups:

  Coordinate/group columns identifying a prediction location.

- value:

  Column containing posterior draw values.

- ...:

  Arguments passed to \`summarise_posterior_matrix()\`.

## Value

A data frame in coordinate order with one row per location.
