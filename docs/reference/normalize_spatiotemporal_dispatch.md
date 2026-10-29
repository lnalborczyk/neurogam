# Normalise the options supplied to the spatiotemporal dispatcher

Normalise the options supplied to the spatiotemporal dispatcher

## Usage

``` r
normalize_spatiotemporal_dispatch(
  sensor_id,
  baseline_window,
  varying_smooth,
  varying_smooth_missing,
  kvalue,
  kvalue_missing,
  bs,
  bs_missing
)
```

## Arguments

- sensor_id:

  Name of the EEG sensor identifier column.

- baseline_window:

  Baseline-ROPE window or NULL.

- varying_smooth:

  Whether participant-specific smooths are requested.

- varying_smooth_missing:

  Whether the public caller omitted that option.

- kvalue, bs:

  Tensor marginal settings supplied by the user.

- kvalue_missing, bs_missing:

  Whether those options were omitted.

## Value

A list containing normalized \`kvalue\` and \`bs\`.
