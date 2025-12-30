<!--

# neurogam 0.0.4

## New features

* Implementing spatio-temporal BGAMs for EEG data.

## Other changes

* ...

## Bug fixes

* ...

-->

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

## Other changes

* Improved documentation for `print()` and `summary()` methods.
* Removed the `multilevel = "full"` option (too slow).
* Now returns clusters with both positive and negative signs in `find_clusters()`.
* Allowing to manually specify fill limits in `plot_eeg()`.
* Allowing `ggplot2` theme to be modified in plotting functions.
* Providing a more informative error message after internal data summary when `outcome_sd` contains NAs.
* Allowing to pass a precomputed outcome SD to `testing_through_time()`.

## Bug fixes

* Fixing aberrant PPCs in the `ppc()` method (and simplifying arguments with `ppc_type = c("group", "participant")`).
* Allowing "negative" clusters in `cluster_results` `print()`, `summary()`, and `plot()` methods.
* Fixing error in `testing_through_time()` when `multilevel = "summary"` and `include_ar_term == TRUE`.

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
