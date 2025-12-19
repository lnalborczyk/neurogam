# Check and visualise residual autocorrelation in a fitted neurogam model

Computes practical diagnostics for residual autocorrelation in a fitted
brms model, with special handling of time-resolved neurogam data where
each participant can have multiple independent time series (e.g., one
per condition).

## Usage

``` r
check_residual_autocorrelation(
  fit,
  data = NULL,
  time_id = "time",
  series_id = "ar_series",
  participant_id = "participant",
  predictor_id = "predictor",
  outcome_mean_id = "outcome_mean",
  success_id = "success",
  trials_id = "trials",
  max_lag = 25,
  n_series_plot = 9,
  seed = 42,
  use_posterior_mean = TRUE,
  verbose = TRUE
)
```

## Arguments

- fit:

  A fitted [`brmsfit`](https://paulbuerkner.com/brms/reference/brm.html)
  object.

- data:

  A data frame used to fit the model (or the same rows in the same
  order). Must contain a time column and (directly or indirectly) an AR
  series identifier. If `NULL` (default), uses data from `fit`.

- time_id:

  Character; name of the time column. Defaults to `"time"`.

- series_id:

  Character; name of the column defining independent time series for
  autocorrelation checks. Defaults to `"ar_series"`. If `series_id` is
  not present in `data`, it is created as
  `interaction(participant_id, predictor_id, drop = TRUE)` when
  possible, otherwise as `participant_id`.

- participant_id:

  Character; participant column name used to build `series` if needed.
  Defaults to `"participant"`.

- predictor_id:

  Character; predictor/condition column name used to build `series` if
  needed. Defaults to `"predictor"`.

- outcome_mean_id:

  Character; outcome column for Gaussian-type residuals. Defaults to
  `"outcome_mean"`.

- success_id, trials_id:

  Character; columns for binomial outcomes (`success` and `trials`).
  Defaults to `"success"` and `"trials"`.

- max_lag:

  Integer; maximum lag for ACF curves. Defaults to 25.

- n_series_plot:

  Integer; number of series to plot ACFs for (sampled). Defaults to 9.

- seed:

  Integer or `NULL`; random seed for sampling series to plot. Defaults
  to 42.

- use_posterior_mean:

  Logical; if `TRUE` (default) uses posterior mean `E[y|x]` to compute
  residuals. If `FALSE`, uses median.

- verbose:

  Logical; if `TRUE`, emits informative messages and warnings.

## Value

A named list with components:

- data:

  A copy of `data` with `.series` and `.resid` columns added.

- checks:

  A list of checks (duplicates, missing columns, etc.).

- ar_summary:

  A data frame summarising AR parameter draws (if present).

- rho1_by_series:

  A data frame with lag-1 residual correlation per series.

- acf_df:

  A data frame of ACF values for sampled series (for plotting).

- plots:

  A list of ggplot2 objects.

## Details

The function:

- checks that time points are unique within each AR series,

- extracts and summarises AR parameter(s) (e.g., `ar[1]`) if present,

- computes posterior-mean residuals \\y - E\[y \mid x\]\\,

- quantifies lag-1 residual autocorrelation per series,

- visualises residual autocorrelation via ACF curves and summary plots.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
# after fitting a brms model
out <- check_residual_autocorrelation(fit)
out$ar_summary
out$plots$rho1_hist
out$plots$acf
} # }
```
