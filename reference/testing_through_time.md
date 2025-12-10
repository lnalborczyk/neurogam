# Time-resolved testing based on BGAMMs

Fits time-resolved Bayesian generalised additive (multilevel) models
(BGAMMs) using brms, and computes posterior odds for an effect at each
time point. The effect can be either i) a deviation of the outcome from
a reference value (e.g., zero or a chance level), or ii) a difference
between two groups/conditions.

## Usage

``` r
testing_through_time(
  data,
  participant_id = "participant",
  outcome_id = "eeg",
  time_id = "time",
  predictor_id = "condition",
  family = gaussian(),
  kvalue = 20,
  bs = "tp",
  multilevel = c("summary", "full", "group"),
  warmup = 1000,
  iter = 2000,
  chains = 4,
  cores = 4,
  backend = "cmdstanr",
  threshold = 10,
  n_post_samples = NULL,
  chance_level = 0,
  sesoi = 0,
  credible_interval = 0.95
)
```

## Arguments

- data:

  A data frame in long format containing time-resolved data.

- participant_id:

  Character; name of the column in `data` specifying participant IDs.

- outcome_id:

  Character; name of the column in `data` containing the outcome values
  (e.g., M/EEG amplitude, decoding accuracy).

- time_id:

  Character; name of the column in `data` containing time information
  (e.g., in seconds or samples).

- predictor_id:

  Character; name of the column in `data` containing a *binary*
  predictor (e.g., group or condition). If `predictor_id = NA`, the
  function tests whether the outcome differs from
  `0 + chance_level + sesoi` over time.

- family:

  A brms family object describing the response distribution to be used
  in the model (defaults to
  [`gaussian()`](https://rdrr.io/r/stats/family.html)).

- kvalue:

  Numeric; basis dimension `k` passed to the smooth term
  `s(time, ..., k = kvalue)`.

- bs:

  Character; Character scalar; type of spline basis to be used by brms
  (passed to `s()`, e.g., `"tp"` for thin-plate splines).

- multilevel:

  Character; which model to fit. One of

  - `"full"`: Full GAMM with participant-level random/varying effects;

  - `"summary"`: GAMM fitted to participant-level summary statistics
    (mean outcome and its standard deviation);

  - `"group"`: Group-level GAM fitted to participant-averaged data (no
    random/varying effects).

- warmup:

  Numeric; number of warm-up iterations per chain.

- iter:

  Numeric; total number of iterations per chain (including warmup).

- chains:

  Numeric; number of MCMCs.

- cores:

  Numeric; number of parallel cores to use.

- backend:

  Character; package to use as the backend for fitting the Stan model.
  One of `"cmdstanr"` (default) or `"rstan"`.

- threshold:

  Numeric; threshold on the posterior odds (`prob_ratio`) used to define
  contiguous temporal clusters. Values greater than 1 favour the
  hypothesis that the effect exceeds `chance_level + sesoi`.

- n_post_samples:

  Numeric; number of posterior draws used to compute posterior
  probabilities. If `NULL` (default), all available draws from the
  fitted model are used.

- chance_level:

  Numeric; reference value for the outcome (e.g., 0.5 for decoding
  accuracy). Only used when testing against a constant (i.e., when there
  is no `predictor_id` or when the effect is a difference from chance).

- sesoi:

  Numeric; smallest effect size of interest (SESOI). The posterior
  probability is computed for the effect being strictly larger than
  `chance_level + sesoi`.

- credible_interval:

  Numeric; width of the credible (quantile) interval.

## Value

An object of class `"clusters_results"`, which is a list with elements:

- `clusters`: a data frame with one row per detected cluster (e.g.,
  `cluster_onset`, `cluster_offset`, `duration`);

- `predictions`: a data frame with time-resolved posterior summaries
  (posterior median, credible interval, posterior probabilities, and
  odds `prob_ratio`);

- `data`: data used to fit the brms model (possibly summarised);

- `model`: the fitted brms model object;

- `multilevel`: the value of the `multilevel` argument.

The object has an associated
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
visualising the smoothed time course and detected clusters.

## Details

Internally, the function:

1.  builds a formula with a smooth term over time (optionally by group);

2.  fits a brms model according to `multilevel`;

3.  uses tidybayes to extract posterior predictions over time;

4.  computes, at each time point, the posterior probability that the
    effect (or condition difference) exceeds `chance_level + sesoi`;

5.  converts this into posterior odds (`prob_ratio`) and applies a
    clustering procedure
    ([`find_clusters()`](https://lnalborczyk.github.io/neurogam/reference/find_clusters.md))
    over time.

## See also

[`brm`](https://paulbuerkner.com/brms/reference/brm.html)

## Author

Ladislas Nalborczyk <ladislas.nalborczyk@cnrs.fr>.

## Examples

``` r
if (FALSE) { # \dontrun{
# import some simulated EEG data
data(eeg_data)
head(eeg_data)

# fit the BGAMM to identify clusters
results <- testing_through_time(data = eeg_data)

# display the identified clusters
print(results$clusters)

# plot the GAM-smoothed signal and identified clusters
plot(results)
} # }
```
