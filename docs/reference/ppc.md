# Posterior predictive checks

Generate posterior predictive checks (PPCs) from a fitted Bayesian
time-resolved GAMM stored in a `clusters_results` object. PPCs can be
produced either at the group level or separately for each participant.

## Usage

``` r
ppc(
  object,
  ppc_type = c("group", "participant"),
  ndraws = 500,
  group_var = NULL,
  cores = 4,
  ...
)
```

## Arguments

- object:

  A `clusters_results` object containing a fitted
  [`brmsfit`](https://paulbuerkner.com/brms/reference/brmsfit-class.html)
  model in `object$model`.

- ppc_type:

  Character string specifying the type of PPC to generate. Either
  `"group"` (default) for group-level PPCs (ignoring participant
  identity) or `"participant"` for participant-wise PPCs.

- ndraws:

  Integer specifying the number of posterior draws to use for the PPC.
  Defaults to 500.

- group_var:

  Optional character; name of the grouping variable to use for grouped
  PPCs at the group level. If NULL (default), the function uses
  "predictor" when present in model\$data and binary (two levels).

- cores:

  Numeric; number of parallel cores to use (only used when
  `ppc_type = "participant"`).

- ...:

  Currently unused. Included for future extensions.

## Value

A `ggplot` object visualising the posterior predictive check. The plot
is printed to the active graphics device and also returned invisibly.

## Details

At the group level, predictions are obtained by simulating from the
posterior using
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
with `re_formula = NA`, after collapsing the original data across
participants (by time). At the participant level, PPCs are generated
using
[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html) with
grouped ribbons.

- **Group-level PPCs** are computed by averaging numeric variables
  across participants at each time point, and simulating posterior
  predictive draws with random effects excluded (`re_formula = NA`).
  This provides a marginal, population-level posterior predictive check.

- **Participant-level PPCs** are computed using grouped ribbon plots,
  showing posterior predictive distributions separately for each
  participant.

The returned object is a `ggplot2` object produced by
[`ppc_ribbon`](https://mc-stan.org/bayesplot/reference/PPC-intervals.html)
or [`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html),
depending on the selected `ppc_type`.

## See also

[`pp_check`](https://mc-stan.org/bayesplot/reference/pp_check.html),
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
[`ppc_ribbon`](https://mc-stan.org/bayesplot/reference/PPC-intervals.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Group-level PPC
ppc(object = res, ppc_type = "group")

# Participant-level PPC
ppc(object = res, ppc_type = "participant")
} # }
```
