# Build a brms formula for time-resolved BGAMMs

Construct the appropriate
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.html)
used internally by
[`testing_through_time`](https://lnalborczyk.github.io/neurogam/reference/testing_through_time.md),
depending on the response distribution (Gaussian vs. Binomial), whether
a categorical predictor is provided, the multilevel mode (`"summary"`
vs. `"group"`), and the chosen random/varying effects structure.

## Usage

``` r
make_bgam_formula(
  family,
  multilevel = c("summary", "group"),
  predictor_type = c("none", "categorical", "continuous"),
  within_between = c(NA, "within-subject", "between-subject"),
  kvalue = 20,
  bs = "tp",
  include_ar_term = FALSE,
  varying_smooth = TRUE,
  include_by_smooth = TRUE
)
```

## Arguments

- family:

  A family object describing the response distribution. Must be either
  [`gaussian()`](https://rdrr.io/r/stats/family.html) or
  [`binomial()`](https://rdrr.io/r/stats/family.html) (as used in
  `testing_through_time`).

- multilevel:

  Character; which model family to build. One of `"summary"`
  (participant-level varying effects) or `"group"` (population-level GAM
  without varying effects).

- predictor_type:

  Character; type of predictor. One of `"none"`, `"categorical"`, or
  `"continuous"`.

- within_between:

  Character; only used when `predictor_type = "categorical"` and
  `multilevel = "summary"`. Must be one of `"within-subject"` or
  `"between-subject"` to determine whether to include a varying slope
  `(1 + predictor | participant)` or only `(1 | participant)`.

- kvalue:

  Numeric; basis dimension `k` used in `s(time, ..., k = kvalue)`.

- bs:

  Character; spline basis used in `s(time, bs = bs, ...)` (e.g.,
  `"tp"`).

- include_ar_term:

  Logical; if `TRUE`, adds an AR(1) autocorrelation structure within
  participant via .
  `autocor = brms::ar(time = "time", gr = "participant", p = 1, cov = FALSE)`.

- varying_smooth:

  Logical; if `TRUE` (default), add a factor-smooth interaction
  `s(participant, time, bs = "fs", ...)` in `multilevel = "summary"`
  models. If `FALSE`, only include the intercept/slope varying effects.

- include_by_smooth:

  Logical; if `TRUE` (default) and `predictor_type = "categorical"`, the
  time smooth is specified with `by = predictor`.

## Value

A
[`brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.html)
object.

## Details

The function assumes that the data passed to `brm()` has already been
reshaped to use the internal column names expected by
`testing_through_time`: `time`, `participant`, `predictor` (optional),
and one of `outcome_mean` / `outcome_sd` (Gaussian) or `success` /
`trials` (Binomial).

## Examples

``` r
if (FALSE) { # \dontrun{
# Gaussian, summary, no predictor
make_bgam_formula(
  family = gaussian(),
  multilevel = "summary",
  predictor_type = "none",
  kvalue = 20, bs = "tp"
  )

# Binomial, summary, within-subject predictor
make_bgam_formula(
  family = binomial(),
  multilevel = "summary",
  predictor_type = "categorical",
  within_between = "within-subject",
  kvalue = 20, bs = "tp"
  )
} # }
```
