# Construct a spatiotemporal Bayesian GAMM formula

Internal formula builder for one temporal and two spatial coordinates.
The spatial coordinates form one marginal of a tensor-product smooth,
and time forms the other marginal.

## Usage

``` r
make_bgam_formula_spatiotemporal(
  family,
  multilevel = "summary",
  predictor_type = "none",
  within_between = NA_character_,
  time_id = "time",
  space_id = c("xproj", "yproj"),
  kvalue = c(15L, 10L),
  bs = c("ps", "ps"),
  tensor_d = c(2L, 1L),
  t2_full = FALSE,
  varying_smooth = FALSE,
  include_ar_term = FALSE,
  include_by_smooth = TRUE,
  use_se = TRUE
)
```

## Arguments

- family:

  A Gaussian family object.

- multilevel:

  Must be \`"summary"\`.

- predictor_type:

  \`"none"\` or \`"categorical"\`.

- within_between:

  \`"within-subject"\` or \`"between-subject"\` when a predictor is
  supplied.

- time_id:

  One temporal variable name.

- space_id:

  Two spatial variable names.

- kvalue:

  Basis dimensions for the two marginal bases.

- bs:

  Basis families for the two marginal bases.

- tensor_d:

  Dimensions per marginal, currently \`c(2, 1)\`.

- t2_full:

  Whether to use the full \`mgcv::t2()\` decomposition.

- varying_smooth:

  Participant-specific 3D smooths are not supported.

- include_ar_term:

  AR residual structure is not supported here.

- include_by_smooth:

  Include a condition-specific smooth.

- use_se:

  Use \`se(outcome_sem, sigma = TRUE)\` on the left-hand side.

## Value

A \`brmsformula\` object.
