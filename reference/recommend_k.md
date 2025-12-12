# Recommend a smooth basis dimension k via effective complexity "knee" detection

This function takes a `"clusters_results"` object produced by
[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
refits the underlying brms model for a grid of `k` values, computes
model comparison criteria (`loo`, `waic`), and identifies a recommended
smooth basis dimension via an automatic "knee" detection procedure on
the effective number of parameters (either `p_loo` or `p_waic`).

## Usage

``` r
recommend_k(
  object,
  k_min = 10,
  k_max = 40,
  k_step = 5,
  criterion = c("waic", "loo"),
  knee_method = c("geometric_smooth", "geometric"),
  loess_span = 0.5,
  verbose = TRUE
)
```

## Arguments

- object:

  An object of class `"clusters_results"` as returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).
  The object must contain a fitted brms model in its `$model` slot.

- k_min:

  Numeric; minimum value of the smooth basis dimension `k` to consider.

- k_max:

  Numeric; maximum value of the smooth basis dimension `k` to consider.

- k_step:

  Numeric; step size between successive `k` values. The sequence of
  tested values is `seq(k_min, k_max, by = k_step)`.

- criterion:

  Character vector of model comparison criteria to add via
  [`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.html).
  Defaults to `c("waic", "loo")`.

- knee_method:

  Character; method for knee detection. One of:

  - `"geometric"`: standard geometric elbow on the raw `p_*` values;

  - `"geometric_smooth"`: geometric elbow on a loess-smoothed curve of
    `p_*` vs `k`, with weights `1 / SE(p_*)^2` when available
    (recommended when `p_*` estimates are noisy).

- loess_span:

  Numeric; smoothing parameter passed to
  [`loess`](https://rdrr.io/r/stats/loess.html) when
  `knee_method = "geometric_smooth"`. Defaults to `0.75`.

- verbose:

  Logical; if `TRUE` (default), prints progress messages while refitting
  models for different `k` values.

## Value

An object of class `"recommend_k_results"`, which is a list with
elements:

- `models`: named list of refitted `brmsfit` objects, one per `k` value;

- `comparison`: data frame summarising model criteria for each `k` (one
  row per `k`), including `p_*` and its SE;

- `k_values`: numeric vector of `k` values that were evaluated;

- `recommended_k`: the `k` value identified as the knee based on the
  chosen effective complexity measure and method;

- `plot`: a `ggplot2` object displaying `p_*` as a function of `k` with
  error bars and the knee highlighted;

- `knee_on`: the complexity measure used (`"p_loo"` or `"p_waic"`);

- `knee_method`: the knee detection method used.

## Details

For each `k` in `seq(k_min, k_max, by = k_step)`, the function:

1.  updates all `s(..., k = ...)` terms in the original brms formula to
    use the new value `k`, and refits the model via
    [`update()`](https://rdrr.io/r/stats/update.html);

2.  calls
    [`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.html)
    to compute the requested criteria, and extracts `p_loo` or `p_waic`
    and their standard errors (where available);

3.  builds a comparison table of criteria values across `k`.

To recommend a basis dimension, the function treats the chosen
complexity measure (`p_loo` or `p_waic`) as a function of `k` and uses
an elbow (knee) method:

- `knee_method = "geometric"` applies the geometric knee procedure
  directly to `p_*` vs `k`;

- `knee_method = "geometric_smooth"` first fits a loess curve `p_* ~ k`
  (with weights `1 / SE(p_*)^2` when available), then applies the
  geometric knee to the smoothed curve. This reduces the influence of
  noisy `p_*` estimates.

The resulting plot shows `p_*` vs `k` with:

- points and a line for the raw `p_*` estimates;

- vertical error bars `± SE(p_*)` (when SEs are available);

- (optionally) a dashed line for the smoothed `p_*` curve;

- a vertical dashed line at the recommended `k` (knee).

## See also

[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
[`brm`](https://paulbuerkner.com/brms/reference/brm.html),
[`add_criterion`](https://paulbuerkner.com/brms/reference/add_criterion.html)
[`loess`](https://rdrr.io/r/stats/loess.html)

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
# import some simulated EEG data
data(eeg_data)
head(eeg_data)

# fit a time-resolved GAMM
res <- testing_through_time(
  data = eeg_data,
  participant_id = "participant",
  outcome_id = "eeg",
  time_id = "time",
  predictor_id = NA,
  kvalue = 20,
  multilevel = "summary"
  )

# recommend an optimal smooth basis dimension k
k_res <- recommend_k(
  object = res,
  k_min = 10,
  k_max = 40,
  k_step = 5,
  criterion = "waic"
  )

# print summary in the console
summary(k_res)

# extract the recommended k value
k_res$recommended_k

# access the comparison table
k_res$comparison

# visualise effective complexity vs. k
k_res$plot
} # }
```
