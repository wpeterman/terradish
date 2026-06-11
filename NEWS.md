terradish 0.0.40 (dev)
---------
* Added a `curvature` argument to `terradish()` and `terradish_algorithm()`
  (`"exact"` default, or `"gauss_newton"`). The Gauss-Newton/Fisher curvature
  drops the two residual-weighted second-derivative terms, is positive
  semidefinite, needs only first derivatives of conductance, and equals the
  exact Hessian at a well-fitting optimum. It flows unchanged into
  `vcov()`/`summary()`/`confint()`, so the Fisher information becomes the
  parameter covariance and the standard errors are information based.
* Added `solver = "block_cg"`, a Jacobi-preconditioned block conjugate gradient
  (O'Leary) sharing one Krylov space across all focal right-hand sides, for
  small or well-conditioned graphs.
* Added `terradish_kron_reduce_tiled()`, an exact tiled (out-of-core) Schur/Kron
  reduction onto the focal sites. Sequential Schur complements compose, so it
  matches `terradish_kron_reduce()` exactly while eliminating the interior tile
  by tile, bounding peak memory to one tile's interior factor plus the sparse
  interface operator.
* Added the DRAGON structured-coalescent directed engine (`R/dragon.R`, Phase 1,
  pure R + Matrix): estimates asymmetric gene flow as a covariate function via a
  directed migration generator and FRAME's Strobeck forward map (expected pairwise
  coalescence times), with an analytic adjoint gradient, three coalescence-rate
  models (uniform, covariate drift surface, FRAME stationary coupling), a
  `dragon()` fit, a `dragon_collinearity()` direction-vs-drift diagnostic, and S3
  methods. The coalescent-correct successor to `terradish_directed()` (which uses a
  direction-blind commute time). Reference-pinned against the validated Python
  prototype (`tests/testthat/test-dragon.R`). Run `devtools::document()` to register
  exports/man pages before `R CMD check`.

terradish 0.0.39
---------
* Clarified Wishart `nu` guidance for microsatellite workflows, AIC/AICc/BIC model-comparison guidance, and CRAN installation text for `corMLPE`
* Added focused regression coverage for non-symbol formula responses and named Hessian output

terradish 0.0.38
---------
* Fixed covariance-response power summaries so no-signal boundary fits count as parameter non-detections, and revised the simulation-design vignette to demonstrate how marker count, site count, and the `tau / sigma` signal ratio affect power

terradish 0.0.37
---------
* Refined settings assessment guidance and related optimizer documentation, expanded simulation-design vignette support files, and updated vignette cross-links while excluding `vignettes/precompute.R` from package builds

terradish 0.0.36
---------
* Removed unsupported and unverifiable formal references from package documentation and vignettes
* Corrected the McCullagh generalized Wishart citation and added formal references for already-cited MLPE and genomic relationship methods

terradish 0.0.35
---------
* Clarified that Wishart `nu` is an effective degrees-of-freedom value, with SNP and microsatellite guidance across help pages, README, and vignettes
* Expanded vignette cross-links for Wishart covariance, spline conductance, and Gaussian scale-of-effect workflows

terradish 0.0.34
---------
* Fixed `covariance_response_power()` so that scenarios in which every replicate failed to converge produce result rows with `NA` conductance correlations instead of aborting the entire cell with "no complete element pairs"; downstream summaries and parameter recovery tables now reflect the failed fits transparently

terradish 0.0.33
---------
* Restored base R visibility behavior for `plot.terradish()`, so top-level plotting calls auto-render while assigned calls remain quiet until explicit `print()`
* Simplified plotting examples in the README and getting-started vignette to match the visible-return plotting contract
* Added focused regression coverage for plotting visibility and explicit print behavior across returned plot objects
* Added `covariance_response_power()` for Wishart covariance-response power screening across focal sample sizes and sampling strategies

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
