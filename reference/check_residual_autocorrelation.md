# Check residual autocorrelation for brms models

Compute and visualise residual autocorrelation diagnostics (ACF and
partial ACF) for `brms` models, optionally comparing raw residuals to
"corrected" residuals obtained by removing the contribution of
autoregressive (AR) terms.

## Usage

``` r
check_residual_autocorrelation(
  fit,
  data = NULL,
  time_id = "time",
  series_id = "ar_series",
  participant_id = "participant",
  predictor_id = "predictor",
  resid_method = "posterior_epred",
  resid_type = c("ordinary", "pearson"),
  ndraws = 1000,
  use_posterior_mean = TRUE,
  ar_summary = c("mean", "median"),
  max_lag = 20,
  n_series_plot = 9,
  seed = 666,
  acf_ndraws = 1000,
  acf_ci = 0.95,
  theme = ggplot2::theme_bw(),
  verbose = TRUE
)
```

## Arguments

- fit:

  A fitted `brmsfit` object.

- data:

  Optional `data.frame` used to compute diagnostics. Defaults to
  `fit$data`. Must match the data used to fit `fit` (same rows and
  order), otherwise residual draws will not align with `data`.

- time_id:

  Name of the time column in `data`. Must be present. The time variable
  is used to order observations within each series prior to computing
  correlations.

- series_id:

  Name of the series (time-series grouping) column in `data`. If
  present, it is used as the grouping variable. If absent, a series
  factor is constructed using `participant_id` and `predictor_id` if
  available (see Details).

- participant_id:

  Name of a participant identifier column. Used to construct a fallback
  series grouping if `series_id` is not present.

- predictor_id:

  Name of a predictor identifier column. Used to construct a fallback
  series grouping if `series_id` is not present.

- resid_method:

  Residual method passed to `residuals.brmsfit()` (argument `method`).
  Defaults to `"posterior_epred"`.

- resid_type:

  Residual type passed to `residuals.brmsfit()` (argument `resid_type`).
  One of `"ordinary"` or `"pearson"`.

- ndraws:

  Number of posterior draws of residuals to request from
  `residuals.brmsfit(..., summary = FALSE)`. These draws are used both
  to compute the pointwise residual summary stored in `data` and to
  compute ACF/PACF bands.

- use_posterior_mean:

  Logical; if `TRUE`, raw residuals stored in `data` are computed as the
  posterior mean across residual draws. If `FALSE`, the posterior median
  is used.

- ar_summary:

  Summary used to form a single vector of AR coefficients for the
  "point" corrected residuals stored in `data`. One of `"mean"` or
  `"median"`.

- max_lag:

  Maximum lag for ACF/PACF computation (positive integer). Lags are
  computed from 1 to `max_lag`.

- n_series_plot:

  Number of series to randomly sample (without replacement) for plotting
  ACF/PACF panels and residual time series panels.

- seed:

  Optional random seed used for sampling series and posterior draws. If
  `NULL`, sampling is not seeded.

- acf_ndraws:

  Number of posterior draws to subsample when computing ACF/PACF
  summaries. If `acf_ndraws > ndraws`, the maximum available draws are
  used.

- acf_ci:

  Credible interval width for ACF/PACF bands (e.g., `0.95` for a 95
  interval). Must lie in `(0, 1)`.

- theme:

  A `ggplot2` theme object added to plots.

- verbose:

  Logical; if `TRUE`, emit messages and warnings about grouping,
  duplicated time values, or posterior draw alignment.

## Value

A list of class `"neurogam_autocor_check"` with components:

- ar_summary:

  A `data.frame` summarising posterior AR parameters (if present).

- rho1_by_series_raw:

  Lag-1 autocorrelation per series computed from `.resid_raw`.

- rho1_by_series_corrected:

  Lag-1 autocorrelation per series computed from `.resid_corrected` (NA
  if no AR parameters).

- acf_raw:

  A `data.frame` containing ACF summaries for raw residuals, with
  columns `.series`, `lag`, `q_low`, `q50`, `q_high`, `point`.

- acf_corrected:

  Same as `acf_raw` but for corrected residuals (empty/NA if no AR).

- pacf_raw:

  A `data.frame` containing PACF summaries for raw residuals, with
  columns `.series`, `lag`, `q_low`, `q50`, `q_high`, `point`.

- pacf_corrected:

  Same as `pacf_raw` but for corrected residuals (empty/NA if no AR).

- plots:

  A named list of `ggplot` objects (ACF/PACF panels, rho1 summaries,
  residual time-series plots, and AR density plots if applicable).

## Details

This function is designed for models fit with `brms` that may include
autoregressive structures (e.g.,
[`ar()`](https://rdrr.io/r/stats/ar.html) terms). Residual
autocorrelation is evaluated per time series (defined by `series_id` or
fallback groupings) and summarised using posterior draws of residuals.

**Series grouping.** If `series_id` is present in `data`, it defines the
grouping. Otherwise, the function tries to construct `.series` as:

1.  `interaction(participant_id, predictor_id)` if both columns exist;

2.  `participant_id` if only participant exists;

3.  `predictor_id` if only predictor exists;

4.  a single-level factor `"group"` if none exist.

**Corrected residuals.** If AR parameters are detected in the posterior
draws (columns matching `"^ar\["`), corrected residuals are computed per
series as: \$\$e_t^{(corr)} = e_t - \sum\_{k=1}^{p} \phi_k e\_{t-k}\$\$
using either the posterior mean or median of \\\phi_k\\ (controlled by
`ar_summary`) for the *point* corrected residuals stored in `data`.

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
# fit is a brmsfit with an AR(1) term
out <- check_residual_autocorrelation(
    fit = fit,
    time_id = "time",
    series_id = "ar_series"
    )

# inspect ACF/PACF summaries
head(out$acf_raw)
head(out$pacf_raw)

# plot ACF/PACF panels
out$plots$acf_raw
out$plots$acf_corrected
} # }
```
