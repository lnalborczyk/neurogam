# Posterior predictive checks for time-resolved BGAMMs

Generates posterior predictive checks (PPCs) for the brms model stored
inside a `"clusters_results"` object, displaying:

1.  a PPC of the predicted outcome distribution (densities), and

2.  a PPC of summary statistics (mean and standard deviation).

## Usage

``` r
ppc_clusters_results(
  object,
  ndraws_time = 100,
  ndraws_stat = 100,
  stat = c("mean", "sd"),
  group_var = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"clusters_results"` as returned by
  [`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md).
  The object must contain a fitted brms model in its `$model` slot.

- ndraws_time:

  Integer; number of posterior draws used for the distribution-level PPC
  (density overlay). Defaults to `100`.

- ndraws_stat:

  Integer; number of posterior draws used for the statistic-level PPC
  (mean/SD). Defaults to `500`.

- stat:

  A character vector of summary statistics to use for the `"stat_2d"`
  PPC. Passed to
  [`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) as
  the `stat` argument. Defaults to `c("mean", "sd")`.

- group_var:

  Optional character scalar specifying the name of the grouping variable
  for the grouped PPC. If `NULL` (default), the function automatically
  uses `"predictor"` when that column exists in the model data and has
  exactly two levels (binary predictor). When a grouping variable is
  used, the density PPC employs `type = "dens_overlay_grouped"`;
  otherwise, `type = "dens_overlay"`.

- ...:

  Currently ignored. Included for future extensibility.

## Value

A patchwork / ggplot2 object containing the combined posterior
predictive checks. The plot is also printed as a side effect.

## Details

This function is a convenience wrapper around
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) and
patchwork. It first tries to detect whether the model contains a binary
grouping variable (by default, a column named `"predictor"` in
`model$data`). If such a variable is found (or if `group_var` is
explicitly provided), the distribution-level PPC is produced using
`type = "dens_overlay_grouped"`, which compares predicted and observed
densities within each group. Otherwise, a standard
`type = "dens_overlay"` PPC is drawn.

The second panel uses `type = "stat_2d"` with the supplied `stat`
argument (by default, mean and standard deviation), allowing inspection
of how well the posterior predictive distribution reproduces key summary
statistics of the data.

The two PPC plots are combined side-by-side using
[`wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html),
and the resulting combined plot is returned (and can be further modified
using standard ggplot2 or patchwork operations).

## See also

[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
[`print.clusters_results`](https://lnalborczyk.github.io/neurogam/reference/print.clusters_results.md),
[`summary.clusters_results`](https://lnalborczyk.github.io/neurogam/reference/summary.clusters_results.md),
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html),
[`wrap_plots`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)

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
  kvalue = 20,
  multilevel = "summary"
  )

# posterior predictive checks (combined plot)
ppc_clusters_results(res)
} # }
```
