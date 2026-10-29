# Fit a BGAM or reuse an existing brms model

Centralizes the brms fitting arguments and the \`previous_model\` path.
This helper does not validate model/data compatibility, construct
prediction grids, or compute any posterior quantities; callers retain
those duties.

## Usage

``` r
fit_or_reuse_bgam(
  formula,
  data,
  family,
  previous_model = NULL,
  n_post_samples = NULL,
  warmup = 1000,
  iter = 2000,
  chains = 4,
  cores = 4,
  threads = NULL,
  backend = "cmdstanr",
  stan_control = NULL,
  file = NULL,
  stan_model_args = NULL,
  extra_args = list(),
  resolve_draws = TRUE,
  reuse_message = "Using `previous_model`: skipping model fitting.",
  fit_fun = brms::brm,
  ndraws_fun = brms::ndraws
)
```

## Arguments

- formula:

  A brms model formula.

- data:

  The model data frame.

- family:

  A brms family.

- previous_model:

  Optional existing \`brmsfit\`.

- n_post_samples:

  Optional requested number of posterior draws.

- warmup, iter, chains, cores, threads:

  Sampling settings.

- backend:

  Backend used by brms.

- stan_control:

  Sampling control options forwarded as \`control\`.

- file:

  Optional model cache filename.

- stan_model_args:

  Optional Stan compilation options. If NULL, uses \`stanc_options =
  list("O1")\` for cmdstanr, and an empty list for rstan.

- extra_args:

  Named or unnamed arguments forwarded to \`brms::brm()\`.

- resolve_draws:

  Whether to populate \`n_post_samples\` from the model if the caller
  supplied NULL. Set FALSE when draw selection is handled separately by
  the caller.

- reuse_message:

  Message printed if an existing fit is reused.

- fit_fun, ndraws_fun:

  Callable dependencies; injectable for fast tests.

## Value

A list containing \`model\`, \`n_post_samples\`, and \`reused\`.
