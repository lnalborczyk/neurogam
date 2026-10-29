# CRAN submission comments

## Submission

This is the first submission of `neurogam` (version 0.1.0) to CRAN.

## Test environments and check results

- macOS Sequoia 15.7.4 (Apple Silicon, aarch64), R 4.5.2: 0 ERRORs, 0 WARNINGs, 1 NOTE.
- Windows Server 2022 (x86_64, UCRT), R-devel (2026-10-28, r90783), via win-builder: 0 ERRORs, 0 WARNINGs, 1 NOTE.
- Ubuntu Linux, via GitHub Actions (CRAN-style R CMD check): checks passed.

Package installation, examples, unit tests, vignette rebuilding, and PDF/HTML manual checks passed on both macOS and Windows R-devel.

## Notes

The incoming feasibility NOTE identifies this as a new submission.

On Windows R-devel, the incoming spelling check also highlighted scientific terms (including 'magnetoencephalography', 'electromyography', 'pupillometry', and 'spatiotemporal') and author surnames; these are intentional terms or names.

The Windows and macOS incoming checks reported an HTTP 429 response from a bioRxiv reference in a precompiled vignette. This is rate limiting by the external website, not a network request made by package code.

## Source package and documentation size

The compressed source archive is approximately 13.5 MB, above the recommended 10 MB. The installed documentation occupies approximately 6.1 MB, above the general 5 MB guideline; the example data occupy approximately 3.9 MB when installed.

The package provides illustrated scientific vignettes demonstrating temporal, temporal-generalisation, and spatiotemporal analysis, together with example electrophysiological datasets. The documentation includes many plots needed to interpret time-resolved results and sensor topographies.

Figure resolution has already been reduced to preserve legibility while minimising size. Bayesian model-fitting computations are precomputed for the distributed vignettes, so these computationally expensive analyses do not run during package checks, and vignette rebuilding is fast.

I would appreciate consideration of this modest size exception because further downsampling or removal of figures would reduce the usefulness of the methodological documentation.
