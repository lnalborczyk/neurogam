# Retrieve chunked posterior expectation draws with paired contrasts

Retrieve chunked posterior expectation draws with paired contrasts

## Usage

``` r
chunked_epred(
  model,
  newdata,
  newdata_b = NULL,
  draw_ids = NULL,
  chunk_size = 256L,
  re_formula = NA,
  predict_fun = brms::posterior_epred,
  summarise_fun = NULL
)
```

## Arguments

- model:

  A fitted model, passed to \`predict_fun\`.

- newdata:

  Prediction grid for condition A.

- newdata_b:

  Optional matched prediction grid for condition B.

- draw_ids:

  Optional indices, identical for every chunk and condition.

- chunk_size:

  Number of prediction locations per chunk.

- re_formula:

  Group-level effects specification.

- predict_fun:

  Function accepting model, newdata, re_formula, draw_ids. Defaults to
  \`brms::posterior_epred\`.

- summarise_fun:

  Optional function applied to each draws-by-location chunk. If NULL,
  return a list of draw matrices.

## Value

List of chunk outputs in prediction-grid order.
