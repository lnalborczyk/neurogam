# neurogam 0.0.4

## New features

* Implementing 2D temporal BGAMs, which can be useful for modelling cross-temporal decoding generalisation matrices. Such models can be fitted by defining two time columns (i.e., training and testing times) in `testing_through_time()`. For instance, `testing_through_time(..., time_id = c(train_time, test_time), ...)`.
* Implementing post-processing methods (e.g., `print()`, `summary()`, `plot()`) for 2D temporal models.
* Adding a `previous_model` argument to `testing_through_time()`. This allows to pass an existing `neurogam` model (previously fitted with `testing_through_time()`), which may be useful for exploring the effects of the `threshold` parameters (this avoids re-fitting the model).

## Other changes

* Including the `timegen_data` (2D temporal data).
* Including a new vignette on 2D temporal models using the `timegen_data` data.
* Including a new vignette on checking and modelling residual auto-correlation (for 1D temporal data).
* Allowing to use within-chain parallelisation in `testing_through_time()` via the `threads` argument.
* Allowing to save the fitted `brms` model in `testing_through_time()` via the `save` argument.
* Re-factoring the `check_residual_autocorrelation()` function to use `residuals.brmsfit()` and returns credible intervals.

## Bug fixes

* Fixing the specification of AR terms in `make_bgam_formula()` (now works for both participant- and group-level models).
* The `check_residual_autocorrelation()` function now works when `multilevel = "group"`.

# neurogam 0.0.3

## New features

* Implementing the `recommend_k()` function, running model comparison with varying basis dimension (k).
* Implementing a `ppc()` method for `cluster_results` objects. (#2)
* Adding support for modelling auto-correlation via `testing_through_time(include_ar_term == TRUE)`.
* Implementing the `check_residual_autocorrelation()` function to visualise auto-correlation of the residuals.
* New `plot_sensors()` and `plot_eeg()` functions to visualise EEG sensors grid, EEG raw data, and GAM predictions.
* Adding support for returning participant-level onsets in `testing_through_time()`, `find_clusters()`, and related methods such as `summary()` and `plot()`. (#4)
* Adding support for `binomial()` responses. (#7)
* Adding support for continuous predictors in `testing_through_time()`. (#1)
* Allowing to pass a precomputed outcome SD to `testing_through_time()`.
* Adding the `use_se` argument in `testing_through_time()` and `make_bgam_formula()`.

## Other changes

* Improved documentation for `print()` and `summary()` methods.
* Removed the `multilevel = "full"` option (too slow).
* Now returns clusters with both positive and negative signs in `find_clusters()`.
* Allowing to manually specify fill limits in `plot_eeg()`.
* Allowing `ggplot2` theme to be modified in plotting functions.
* Providing a more informative error message after internal data summary when `outcome_sd` contains NAs.

## Bug fixes

* Fixing aberrant PPCs in the `ppc()` method (and simplifying arguments with `ppc_type = c("group", "participant")`).
* Allowing "negative" clusters in `cluster_results` `print()`, `summary()`, and `plot()` methods.
* Fixing error in `testing_through_time()` when `multilevel = "summary"` and `include_ar_term == TRUE`.
* Fixing the `varying_smooth_term` error in `make_bgam_formula()` when `predictor_id = NA`.

# neurogam 0.0.2

## New features

* Allowing 3 different models to be fitted: full GAMM, GAMM with summary statistics (recommended), or group-level GAM.
* Adding further support for presence or absence of predictor (e.g., group, condition). When `predictor_id = NA`, `neurogam` now tests signal against 0 through time.
* Implementing `print()` and `summary()` methods for `cluster_results` objects.
* Improved plotting: now plotting the GAM predictions with raw data and improved clusters aesthetics.

## Other changes

* Improved functions documentation and new package website.
* Factoring posterior odds computation within internal functions.

## Bug fixes

* Fixing group-level posterior predictions when multilevel is "full" or "summary".

# neurogam 0.0.1

* Pushing the first version of `neurogam`.
