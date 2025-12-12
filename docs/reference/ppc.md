# Posterior predictive checks for `clusters_results` objects

Generates a posterior predictive check for the brms model stored inside
a `"clusters_results"` object. The function is a thin wrapper around
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) that
automatically chooses a grouping variable when possible.

## Usage

``` r
ppc(
  object,
  prefix = c("ppc", "ppd"),
  ppc_type = "ribbon_grouped",
  ndraws = 100,
  group_var = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"clusters_results"` as returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).
  The object must contain a fitted brms model in its `$model` slot.

- prefix:

  Character; passed to
  [`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) to
  control whether the function uses posterior predictive checks
  (`"ppc"`) or posterior predictive distributions (`"ppd"`). One of
  `"ppc"` (default) or `"ppd"`.

- ppc_type:

  Character; the type of check passed to
  [`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  (default: `"ribbon_grouped"`). See brms documentation for available
  PPC/PPD types.

- ndraws:

  Numeric; number of posterior draws used to generate the PPC/PPD
  (default: `100`).

- group_var:

  Optional character; name of the grouping variable to use for grouped
  PPCs. If `NULL` (default), the function uses `"predictor"` when
  present in `model$data` and binary (two levels). Otherwise it falls
  back to `"participant"`.

- ...:

  Additional arguments passed to
  [`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html).

## Value

A single ggplot2 object (invisibly). The plot is also printed as a side
effect.

## Details

If `group_var` is `NULL`, the function attempts to detect whether the
model data contain a binary grouping variable named `"predictor"`. If
so, it produces a grouped PPC using `group = "predictor"`; otherwise it
uses `group = "participant"`.

The check is performed over `x = "time"` and uses `re_formula = NA` when
the grouping variable is `"predictor"` (i.e., group-level PPC). When
falling back to participant-level grouping, random effects are included
by default (unless overridden via `...`).

## See also

[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# import some simulated EEG data
data(eeg_data)
head(eeg_data)

# fit time-resolved model (one-sample test)
res <- testing_through_time(
  data = eeg_data,
  participant_id = "participant",
  outcome_id = "eeg",
  time_id = "time",
  predictor_id = NA,
  kvalue = 10,
  multilevel = "summary"
  )

# posterior predictive checks
ppc(res)
} # }
```
