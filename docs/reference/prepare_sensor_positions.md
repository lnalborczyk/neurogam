# Prepare a sensor coordinate table

Prepare a sensor coordinate table

## Usage

``` r
prepare_sensor_positions(
  data,
  sensor_id = "sensor",
  spatial_id = c("xproj", "yproj")
)
```

## Arguments

- data:

  Data frame with one or more rows per sensor.

- sensor_id:

  Name of the sensor identifier column.

- spatial_id:

  Names of two numeric projected spatial columns.

## Value

One row per sensor, sorted by sensor identifier.
