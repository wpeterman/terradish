terradish 0.0.25
---------
* Added covariance-response marginal plots for `wishart_covariance` fits, including support for `smooth_loglinear_conductance`
* Updated plot defaults with cleaner axis styling, compact numeric labels, and cowplot-inspired panel presentation without adding a cowplot dependency
* Added `marginal_covariates` to subset marginal and marginal-response plot panels
* Updated cross-validation helpers to accept supported conductance model factories instead of hardcoding `loglinear_conductance`
* Rebuilt fixed-graph conductance factories on train, test, and full cross-validation surfaces when supported
* Expanded focused tests for plotting, Wishart covariance responses, smooth conductance, Gaussian scale conductance, and cross-validation helpers

terradish 0.0.24
---------
* Added spline-based log-linear conductance workflow with `smooth_loglinear_conductance`
* Added a tutorial vignette for spline conductance models
* Added plotting metadata so spline marginal effects are shown against original covariate names

terradish 0.0.23
---------
* Improved optimizer resilience by switching failed Hager-Zhang line searches to bounded backtracking
* Added diagnostics for line-search fallback behavior

terradish 0.0.22
---------
* Added modal imputation for missing allele calls in genetic covariance utilities
* Preserved imputation details as attributes for downstream inspection

terradish 0.0.21
---------
* Stopped tracking local future-development notes in package history
* Kept local planning notes separate from package source files

terradish 0.0.20
---------
* Generalized genetic covariance utilities for individual- and population-level workflows
* Added helpers for converting covariance matrices to pairwise distances
* Improved documentation for covariance-response modeling workflows

terradish 0.0.19
---------
* Added selected-pair measurement models for fitting subsets of pairwise responses
* Added support infrastructure for endpoint covariates in pairwise measurement models

terradish 0.0.18
---------
* Hardened Windows installation workflows
* Refreshed README guidance for package setup and usage

terradish 0.0.17
---------
* Added optimization setting assessment tools
* Added summaries for solver and approximation choices

terradish 0.0.16
---------
* Added cached CHOLMOD direct backend support
* Improved reuse of sparse direct solver setup across repeated evaluations

terradish 0.0.15
---------
* Added large-raster efficiency workflows
* Added approximation support for larger landscape genetic surfaces

terradish 0.0.14
---------
* Improved CRAN namespace compliance
* Cleaned package exports and documentation metadata

terradish 0.0.13
---------
* Added Gaussian scale plotting support
* Refreshed package documentation for Gaussian scale-aware conductance workflows
