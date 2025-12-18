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
  trials_id = NULL,
  family = gaussian(),
  kvalue = 20,
  bs = "tp",
  multilevel = c("summary", "group"),
  participant_clusters = FALSE,
  varying_smooth = TRUE,
  warmup = 1000,
  iter = 2000,
  chains = 4,
  cores = 4,
  backend = c("cmdstanr", "rstan"),
  stan_control = NULL,
  n_post_samples = NULL,
  threshold = 10,
  threshold_type = c("both", "above", "below"),
  chance_level = NULL,
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

  Character; name of the column in `data` containing either:

  - A *binary* categorical predictor (e.g., group or condition), in
    which case the function tests, at each time point, whether the
    difference between the two levels exceeds `chance_level`;

  - A *continuous* numeric predictor, in which case the function tests,
    at each time point, whether the *slope* of the outcome with respect
    to the predictor differs from `chance_level` (typically with
    `chance_level = 0`).

  - If `predictor_id = NA`, the function tests whether the outcome
    differs from `chance_level` over time (useful for decoding
    accuracies, for instance).

- trials_id:

  Character; name of the column in `data` containing the number of
  trials when using `family = binomial()` and summary data. If NULL
  (default), the function internally summarise binary data into
  "successes" and total number of "trials".

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

  - `"summary"`: GAMM fitted to participant-level summary statistics
    (mean outcome and its standard deviation);

  - `"group"`: Group-level GAM fitted to participant-averaged data (no
    random/varying effects).

- participant_clusters:

  Logical; should we return clusters at the participant-level.

- varying_smooth:

  Logical; should we include a varying smooth. Default is TRUE. If
  FALSE, we only include a varying intercept and slope.

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

- stan_control:

  List; parameters to control the MCMC behaviour, using default
  parameters when NULL. See `?brm` for more details.

- n_post_samples:

  Numeric; number of posterior draws used to compute posterior
  probabilities. If `NULL` (default), all available draws from the
  fitted model are used.

- threshold:

  Numeric; threshold on the posterior odds used to define contiguous
  temporal clusters. Values greater than 1 favour the hypothesis that
  the effect exceeds `chance_level`.

- threshold_type:

  Character scalar controlling which clusters are detected. Must be one
  of `"above"`, `"below"`, or `"both"` (default). When `"above"`,
  clusters are formed where `value >= threshold`. When `"below"`,
  clusters are formed where `value <= 1/threshold`. When `"both"`, both
  types are detected and the returned data include a `sign` column.

- chance_level:

  Numeric; null value for the outcome (e.g., 0.5 for decoding accuracy).

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
visualising the smoothed time course and detected clusters, as well as
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods.

## Details

Internally, the function:

1.  builds a formula with a smooth term over time (optionally by group);

2.  fits a brms model according to `multilevel`;

3.  uses tidybayes to extract posterior predictions over time;

4.  computes, at each time point, the posterior probability that the
    effect (or condition difference) exceeds `chance_level`;

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
summary(results)

# plot the model predictions and identified clusters
plot(results)
} # }
```
