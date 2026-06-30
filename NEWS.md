terradish 0.0.44 (dev)
---------
* Fixed the Laplace marginal log-likelihood in `terradish_hierarchical()` to
  use the joint penalized curvature of the covariate and field coefficients.
  The previous field-block-only determinant could overstate support for large
  `tau2` values when the field started absorbing covariate signal.
* Fixed an infinite loop in the Hager-Zhang line search (`bisect()` in
  `hager_zhang.R`), used by `BoxConstrainedNewton()` for every measurement model's
  nuisance-parameter fit. The bisection terminated on an absolute interval-width
  threshold (`.Machine$double.eps`), but the smallest representable gap between two
  doubles of magnitude `|alpha|` is `~|alpha| * eps`; once the bracketing step
  lengths grew past 1, the midpoint `(a + b) / 2` could no longer fall strictly
  between `a` and `b`, the interval stalled above the absolute threshold, and the
  loop spun forever. A flat, numerically noisy nuisance objective (seen with
  `generalized_wishart` on some fits) pushes the line search into exactly that
  regime. Now uses a scale-relative threshold (`eps * max(1, |a|, |b|)`) plus a
  midpoint-stagnation guard; small-step behavior is unchanged.

terradish 0.0.43 (dev)
---------
* Fixed `terradish_directed_algorithm()` to normalize the directed commute-time
  covariance `E` before the measurement-model subproblem (dividing `dL/dE` back by
  the same factor). On large graphs `E` grows with the graph (commute time is
  roughly `2 * |edges| *` resistance), which ill-conditioned the nuisance-parameter
  Hessian and made `BoxConstrainedNewton()`'s diagonal-eigenvalue inverse fail with
  "system is computationally singular"; the optimizer's `tryCatch` swallowed the
  error and the objective collapsed to the `1e12` sentinel, so directed fits
  silently failed to move off the start point on fine rasters. Because the
  measurement-model likelihood enters `E` linearly, the rescaling is exactly
  likelihood-invariant: the log-likelihood, the directional `gamma`, and the
  gradient are unchanged (verified identical across normalization scales to
  ~1e-18). Directed fits now run out-of-the-box at landscape scale with the
  reference `"matrix"` solver.

terradish 0.0.42 (dev)
---------
* Added `solver = "sparse_lu_cpp"` for `terradish_directed()` and
  `terradish_directed_algorithm()`, using an Eigen SparseLU C++ backend that
  caches one factorization per focal absorber and reuses it for the forward
  hitting-time solve and transpose adjoint solve.
* Added regression coverage comparing directed SparseLU covariance and gradient
  output against the reference Matrix backend, and updated the directional
  vignette performance guidance.
* Added the lowercase `logml` return field to `terradish_hierarchical()` as an
  alias for the existing `logML` marginal-likelihood value, including fixed
  numeric `tau2` fits, and clarified that `loglik` excludes the field penalty.

terradish 0.0.41 (dev)
---------
* Added `directed_rates()` and new `plot.terradish_directed()` directional
  visualizations for edge-rate bias and combined symmetric-conductance plus
  directional-bias maps.
* Expanded the directional-conductance vignette with a cached fitted example,
  likelihood-ratio comparison against the reversible special case, edge-rate
  summaries, and interpretation guidance for directional plots.
* Added a focused large-landscape curvature benchmark script and corresponding
  vignette guidance comparing exact and Gauss-Newton curvature.

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
* Added `terradish_kron_reduce_tiled()`, an exact out-of-core Schur/Kron
  reduction onto the focal sites that matches `terradish_kron_reduce()` exactly
  while keeping peak memory in check. By default (`method = "auto"`) it uses the
  fast single-shot reduction while its estimated factorization fits `mem_budget`
  (default 4 GB), and falls back to recursive nested dissection only when it
  would not -- so the common case stays fast and the bounded path is reserved
  for the regime where the single-shot factor (and a one-level tiling, whose
  separator factorization grows ~N^1.6) would run out of memory. Nested
  dissection bisects the interior by a thin separator into two independent
  halves, reduces each recursively, and eliminates the separator last (a
  multifrontal extend-add), bounding every single factorization by a separator
  width; `n_tiles` is its leaf size. With `cores > 1` the independent recursive
  halves are reduced in parallel by forking on Unix (the budget halves down the
  recursion); the chosen `method` and the factor estimate are returned. Supplying
  an explicit `tiles` partition instead uses a flat two-level substructuring that
  parallelizes across tiles (fork on Unix, socket cluster on Windows). Every path
  is identical to the sequential single-shot result, and focal
  vertices are never eliminated, so focal effective resistances are preserved
  exactly.
* Extracted the shared landscape-genetic primitives to the new `landgraph` package
  and now import them from there: the genetic covariance/distance helpers
  (`cov_from_biallelic()`, `cov_from_genetic_data()`, `fst_from_biallelic()`,
  `dist_from_cov()`, `dist_from_biallelic()`) and the directional edge-covariate
  builders (`edge_gradient()`, `edge_flow()`). These are re-exported from terradish,
  so existing code and documentation links are unchanged.
* Moved the DRAGON structured-coalescent directed engine to its own package,
  `dragonflow` (asymmetric gene flow as a function of directional covariates). It is
  no longer shipped in terradish; install `dragonflow` to use `dragon()` and
  `dragon_collinearity()`. terradish remains focused on symmetric resistance.

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
