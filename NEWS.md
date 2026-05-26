terradish 0.0.32
---------
* Clarified that `pair_subset_measurement_model()` does not support full-matrix Wishart likelihoods
* Pointed `generalized_wishart` and `wishart_covariance` users to `wishart_covariates()` for site-level covariance-kernel workflows

terradish 0.0.31
---------
* Expanded vignette guidance for Wishart kernel covariates in joint IBR and IBE workflows
* Stabilized generalized Wishart kernel-covariate fits by symmetrizing the projected inverse before eigendecomposition
* Added a regression test for finite, real-valued generalized Wishart kernel-covariate fits

terradish 0.0.30
---------
* Added `wishart_covariates()` for Wishart measurement models with nonnegative site-level covariance-kernel weights
* Supported both `wishart_covariance` and `generalized_wishart` likelihoods through the new kernel-covariate measurement-model factory
* Added focused tests for kernel construction, base Wishart equivalence, generalized Wishart support, subsetting, and optimizer integration

terradish 0.0.29
---------
* Added the package sticker to the main terradish help pages
* Moved legacy `radish*` wrapper documentation to internal help topics while preserving compatibility methods and aliases
* Added a regression test confirming `summary()` still prints compact terradish summaries when legacy `radish` classes are absent

terradish 0.0.28
---------
* Updated vignettes so plots are explicitly printed under the silent plotting default
* Fixed the IBE/IBR vignette to use a stable pairwise altitude covariate name
* Render-checked all package vignettes to confirm plotting output appears as expected

terradish 0.0.27
---------
* Added package-native focal-support clamping controls for plotting and conductance prediction via `support`, `support_probs`, and `clamp_covariates`
* Extended support-constrained marginal plotting to Gaussian scale-aware conductance workflows
* Updated plotting defaults to return per-panel outputs as individual plots (single `ggplot` or named list) and keep plotting calls silent unless explicitly printed
* Expanded regression tests for support clamping in standard and Gaussian marginal workflows and support-constrained conductance prediction
* Updated README and core vignettes with documented support-clamping usage patterns for stable tail interpretation

terradish 0.0.26
---------
* Added conditional confidence intervals for nuisance parameters in `summary()`
* Added `mlpe_response_change()` to summarize MLPE pairwise-covariate effects on the genetic-distance response scale
* Documented and demonstrated response-scale MLPE covariate interpretation in the IBE/IBR workflow

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
