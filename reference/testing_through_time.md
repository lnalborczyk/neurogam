# Time-resolved testing based on BGAMMs

Fits time-resolved Bayesian generalised additive (multilevel) models
(BGAMMs) using brms, and computes posterior odds for an effect at each
time point. The effect can be either i) a deviation of the outcome from
a reference value (e.g., zero or a chance level), ii) a difference
between two groups/conditions (varying within or between participants),
or iii) the amplitude of a continuous predictor varying either within
(e.g., speech formants) or between participants (e.g., age).

## Usage

``` r
testing_through_time(
  data,
  previous_model = NULL,
  participant_id = "participant",
  outcome_id = "eeg",
  outcome_sem = NULL,
  time_id = "time",
  predictor_id = "condition",
  trials_id = NULL,
  family = gaussian(),
  kvalue = 20,
  bs = "tp",
  multilevel = c("summary", "group"),
  include_ar_term = FALSE,
  use_se = TRUE,
  t2_full = FALSE,
  participant_clusters = FALSE,
  varying_smooth = TRUE,
  include_by_smooth = TRUE,
  include_by_varying_smooth = TRUE,
  warmup = 1000,
  iter = 2000,
  chains = 4,
  cores = 4,
  threads = NULL,
  backend = getOption("brms.backend", "rstan"),
  stan_control = NULL,
  file = NULL,
  n_post_samples = NULL,
  threshold = 10,
  threshold_type = c("both", "above", "below"),
  odds_mode = c("directional_rope", "legacy_one_sided", "legacy_directional",
    "spatiotemporal"),
  chance_level = NULL,
  credible_interval = 0.95,
  baseline_window = c(0, 0.1),
  space_id = NULL,
  sensor_id = NULL,
  tensor_d = c(2L, 1L),
  sensor_adjacency = NULL,
  adjacency_k = 4L,
  adjacency_distance_factor = 2,
  rope = NULL,
  contrast_levels = NULL,
  chunk_size = 256L,
  re_formula = NA,
  ...
)
```

## Arguments

- data:

  A data frame in long format containing time-resolved data.

- previous_model:

  Optional. A previously fitted `brmsfit` object obtained from
  `testing_through_time()`. When provided, the model is not refitted;
  instead, posterior predictions and inference are recomputed using the
  supplied model. This is useful for exploring the effect of different
  `threshold` and `threshold_type` values without re-running model
  fitting.

  The supplied model must be compatible with the current function call
  (i.e., same data structure, formula, family, and predictors). If
  `previous_model` is not `NULL`, arguments related to model estimation
  (e.g., `warmup`, `iter`, `chains`, `cores`, `backend`, `stan_control`)
  are ignored.

- participant_id:

  Character; name of the column in `data` specifying participant IDs.

- outcome_id:

  Character; name of the column in `data` containing the outcome values
  (e.g., M/EEG amplitude, decoding accuracy).

- outcome_sem:

  Character; name of the column in `data` containing the outcome SEM,
  when `outcome_id` has already been summarised (default value is NULL).

- time_id:

  Character; name of the column(s) in `data` containing time information
  (e.g., in seconds or samples).

- predictor_id:

  Character; name of the column in `data` containing either:

  - A *binary* categorical predictor (e.g., group or condition), in
    which case the function tests, at each time point, whether the
    difference between the two levels differs from `chance_level`;

  - A *continuous* numeric predictor, in which case the function tests,
    at each time point, whether the difference between the average value
    of the predictor +1 SD and the average value -1 SD differs from
    `chance_level` (typically with `chance_level = 0`).

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

- include_ar_term:

  Logical; if `TRUE`, adds an AR(1) autocorrelation structure within
  participant via
  `autocor = brms::ar(time = "time", gr = "participant", p = 1, cov = FALSE)`.

- use_se:

  Logical; whether to include known or internally computed measurement
  error via `y | se(outcome_sem)` in the model formula.

- t2_full:

  Logical; If TRUE, then there is a separate penalty for each
  combination of null space column and range space, see
  [`t2`](https://rdrr.io/pkg/mgcv/man/t2.html). Only use when fitting 2D
  temporal models (i.e., when `time_id` contains two temporal
  variables).

- participant_clusters:

  Logical; should we return clusters at the participant-level.

- varying_smooth:

  Logical; should we include a varying smooth. Default is `TRUE`. If
  `FALSE`, we only include a varying intercept and slope.

- include_by_smooth:

  Logical; if `TRUE` (default) and `predictor_type = "categorical"`, the
  time smooth is specified with `by = predictor`.

- include_by_varying_smooth:

  Logical; if `TRUE` (non default) and `predictor_type = "categorical"`,
  the time varying smooth is specified with `by = predictor`.

- warmup:

  Numeric; number of warm-up iterations per chain.

- iter:

  Numeric; total number of iterations per chain (including warmup).

- chains:

  Numeric; number of MCMCs.

- cores:

  Numeric; number of parallel cores to use.

- threads:

  Numeric; number of threads to use in within-chain parallelisation. See
  [`brm`](https://paulbuerkner.com/brms/reference/brm.html)
  documentation for more information.

- backend:

  Character; package to use as the backend for fitting the `Stan` model.
  One of `"cmdstanr"` (default) or `"rstan"`.

- stan_control:

  List; parameters to control the MCMC behaviour, using default
  parameters when NULL. See `?brm` for more details.

- file:

  Either NULL or a character string. In the latter case, the fitted
  `brms` model object is saved via saveRDS in a file named after the
  string supplied in file. The `.rds` extension is added automatically.
  If the file already exists, `brm` will load and return the saved model
  object instead of refitting the model.

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

- odds_mode:

  Posterior-odds inference strategy. The default,
  \`"directional_rope"\`, computes independent continuous posterior odds
  for effects above and below the ROPE. These odds are thresholded
  separately to identify positive and negative clusters.
  \`"legacy_one_sided"\` reproduces the historical capped one-sided
  odds. \`"legacy_directional"\` and \`"spatiotemporal"\` are retained
  for backward compatibility.

- chance_level:

  Numeric; null value for the outcome (e.g., 0.5 for decoding accuracy).

- credible_interval:

  Numeric; width of the credible (quantile) interval.

- baseline_window:

  Numeric; bounds of the baseline window (in the same units of
  `time_id`), used to compute the baseline ROPE. Typically, this is a
  100-200 ms window preceding stimulus onset.

- space_id:

  Optional length-two spatial coordinate names. Non-NULL selects the
  spatiotemporal Gaussian multilevel workflow.

- sensor_id:

  Sensor column name, required when \`space_id\` is supplied.

- tensor_d:

  Tensor marginal dimensions; currently \`c(2, 1)\`.

- sensor_adjacency:

  Optional named sensor adjacency matrix or edge list.

- adjacency_k, adjacency_distance_factor:

  Parameters for automatic sensor adjacency construction when no graph
  is supplied.

- rope:

  Optional nonnegative symmetric ROPE half-width. The default (\`NULL\`)
  estimates the ROPE from the maximum absolute posterior-median
  deviation within \`baseline_window\` for 1D temporal analyses. An
  explicit value, including zero, overrides baseline estimation. For 2D
  and spatiotemporal analyses, \`NULL\` currently corresponds to a
  point-null ROPE of zero.

- contrast_levels:

  Two condition labels specifying the A - B contrast.

- chunk_size:

  Number of sensor/time rows per posterior prediction chunk.

- re_formula:

  Group-level effects in posterior predictions; by default \`NA\`
  excludes subject-specific random effects.

- ...:

  Additional arguments passed to
  [`brm`](https://paulbuerkner.com/brms/reference/brm.html), and from
  there to the selected Stan backend. This can be used, for example, to
  set `refresh`, `seed`, `init`, or other backend-specific sampling
  arguments. Arguments controlling the model formula, data, family,
  backend, MCMC settings, and Stan compilation are handled explicitly by
  `testing_through_time()` and cannot be passed through `...`.

## Value

An object of class `"clusters_results"`, which is a list with elements:

- `clusters`: a data frame with one row per detected cluster (e.g.,
  `cluster_onset`, `cluster_offset`, `duration`);

- `predictions`: a data frame with time-resolved posterior summaries
  (posterior median, credible interval, posterior probabilities, and
  odds `prob_ratio`);

- `model`: the fitted brms model object;

- `summary_data`: data used to fit the brms model (possibly summarised);

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

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a small, reproducible EEG-like example
set.seed(123)

eeg_data <- expand.grid(
    participant = sprintf("participant_%02d", seq_len(10) ),
    condition = c("cond1", "cond2"),
    trial = seq_len(length.out = 5),
    time = seq(from = -0.2, to = 0.6, by = 0.01),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
    )

eeg_data <- eeg_data |>
    mutate(
        eeg = 0.8 * exp(-((time - 0.3) / 0.09)^2) * (condition == "cond2") +
            rnorm(nrow(eeg_data), sd = 0.5)
    )

# fit the BGAMM to identify clusters
results <- testing_through_time(data = eeg_data)

# display the identified clusters
summary(results)

# plot the model predictions and identified clusters
plot(results)

# posterior predictive checks
ppc(results)
} # }
```
