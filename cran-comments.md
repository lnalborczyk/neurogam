# CRAN submission comments

## Resubmission

This is a resubmission of `neurogam` (version 0.1.1), following
CRAN feedback on the initial submission (version 0.1.0).

Changes made in response to the review:

- Removed the redundant "Provides tools for" wording from
  the DESCRIPTION.
- Revised documentation examples, unwrapping fast examples.
- Retained `\dontrun{}` for computationally intensive Bayesian
  model-fitting examples, as agreed with the CRAN reviewer,
  and added explanatory comments to these examples.
- Updated the package to version 0.1.1, incorporating minor
  improvements since the initial submission.

## Test environments and check results

- macOS Sequoia 15.7.4 (Apple Silicon, aarch64), R 4.5.2: 0 ERRORs, 0 WARNINGs, 0 NOTEs.
- Ubuntu Linux (x86_64), R release, GitHub Actions full CRAN-style check (--as-cran --no-manual): 0 ERRORs, 0 WARNINGs; check completed successfully.

The package was checked locally using `devtools::check(cran = TRUE)`.

Package installation, documentation, examples, unit tests, vignettes, and vignette rebuilding all passed.

## Source package and documentation size

The compressed source archive was approximately 13.5 MB in the previous submission, exceeding the recommended 10 MB limit. The latest check reports an installed package size of 6.8 MB, including 5.5 MB of documentation.

The package provides illustrated scientific vignettes demonstrating temporal, temporal-generalisation, and spatiotemporal analyses, together with example electrophysiological datasets. The documentation includes numerous figures needed to interpret time-resolved results and sensor topographies.

Figure resolution has already been reduced to preserve legibility while minimising size. Bayesian model-fitting computations are precomputed for the distributed vignettes, avoiding expensive analyses during package checks and allowing fast vignette rebuilding.

I would appreciate consideration of this modest size exception, as further downsampling or removal of figures would reduce the usefulness of the methodological documentation.
