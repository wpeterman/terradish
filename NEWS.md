terradish 0.0.46
---------
* Reconciled the package overview, function help, examples, vignettes, and
  benchmark guidance with the shared process-model and measurement-model
  hierarchy. Documentation now distinguishes conditional conductance
  associations from causal movement or demographic interpretations, states the
  admissibility requirements for generalized-Wishart distances, separates the
  roles of Wishart `tau`, `sigma`, and `nu`, and treats Gaussian `sigma` as a
  raster-smoothing parameter rather than a movement scale.
* Added comparison-contract checks to `aic_table()` and `anova()` so models
  with different responses, likelihood families, pair subsets, or Wishart
  degrees of freedom are rejected. `AIC.terradish()` now honors its `k`
  argument, and information-criterion tables retain the correct log likelihood
  when model labels are duplicated.
* `simulate_covariance_response()` now rejects Wishart degrees of freedom below
  the covariance dimension and preserves valid noninteger values instead of
  silently truncating them.
* Replaced examples that applied `generalized_wishart` directly to an unchecked
  FST matrix with admissible simulated distance responses or Gaussian-distance
  likelihoods.
* Added `check_distance_response()` and automatic fit-time validation for
  generalized-Wishart responses. The diagnostic tests the positive
  semidefiniteness of the centered Gram matrix, reports the magnitude of any
  violation, and rejects inadmissible or malformed matrices without silently
  applying a Euclidean correction.
* Trimmed the vendored `amgcl` header library from 137 files (1.56 MB) to the
  25 headers `src/` actually compiles against (312 KB), removing the MPI, CUDA,
  VexCL, Eigen, and Epetra backends the package never builds. The kept set is
  the transitive `#include` closure of the seven headers included from
  `src/*.cpp`; the package compiles and the AMG solver tests pass unchanged.
  This also relieves a Windows `MAX_PATH` hazard, since the deleted
  `inst/include/amgcl/mpi/...` paths were the deepest in the tarball.
* `configure` and `cleanup` are now tracked with the execute bit set, so a Unix
  build straight from a git clone works without `R CMD build` having to correct
  the mode.
* Replaced the remaining `1:length()`, `1:nrow()`, and `1:maxit` loop bounds
  with `seq_along()`, `seq_len()`, and `sample.int()`. The optimizers now
  validate `maxit` up front instead of relying on `1:maxit` to produce a
  sensible sequence, so `maxit = 0` errors clearly rather than silently
  iterating twice.
* Documented the last two undocumented S3 methods,
  `print.terradish_covariance_power` and `print.terradish_setting_assessment`,
  each with a "how to read the output" walkthrough, `\seealso`, and a runnable
  example. All 47 registered S3 methods now have a help page. Added examples
  and `\seealso` to `?terradish_methods` and
  `?terradish_cv_replicates_methods`, which previously had neither; the latter
  also gained guidance on reading the mean and standard deviation of the
  held-out log-likelihood.
* `terradish()` gains a `verbose` argument and **fitting is now quiet by
  default**. Previously the default `control` carried `verbose = TRUE`, so every
  fit printed a Newton-Raphson iteration trace, and it did so with `cat()`, which
  `suppressMessages()` cannot silence. All 27 optimizer and line-search traces
  now go through `message()`. `verbose = TRUE` restores the per-iteration
  report; passing it explicitly overrides `control$verbose`, while leaving it
  alone lets a hand-built `control` object speak for itself. The noisier
  line-search trace stays a separate opt-in through
  `NewtonRaphsonControl(ls.control = HagerZhangControl(verbose = TRUE))`.
  `terradish_hierarchical()`, `terradish_assess_settings()`, and
  `terradish_scale_optim()` also default to `verbose = FALSE` now.
* `terradish_hierarchical()` no longer aborts when one `tau2_grid` value fails.
  A large tau2 weakens the field penalty and can leave the nuisance subproblem
  numerically singular, so following the package's own "consider widening
  `tau2_grid`" warning could kill the fit outright. Failing grid points are now
  skipped with a warning naming them, recorded as `NA` in the `tau2_selection`
  table, and excluded from the selection; the grid-edge warning now describes
  the edge of the *usable* range. A fixed `tau2` that fails reports what
  happened and what to try instead, rather than surfacing a raw LAPACK error.
* Every `\dontrun{}` example is gone. The eleven that remained referenced
  objects that were never defined, or were self-contained but slow, so none of
  them could be copied and run. They are now complete, runnable `\donttest{}`
  examples on the bundled `melip` data or a small synthetic lattice, each
  verified to run in a few seconds (the whole set takes about 65 seconds).
  `terradish_scale_optim()` also gained a much fuller `@return`, and its
  `lower`/`upper` documentation now states that the default `scale_fun` measures
  sigma in raster **cells**, not map units; the old example's bounds were far
  outside a useful range.

* Fixed the profile-likelihood Hessian when a nuisance parameter is estimated on
  an active box constraint. `radish_subproblem()` applied the implicit-function
  correction `dphi/dE` to every nuisance parameter, including ones pinned at a
  bound, which cannot move. The correction is now restricted to the free
  parameters. This affected `wishart_covariates()` whenever a kernel
  coefficient was estimated at exactly 0: the reported Hessian, and therefore
  every standard error, confidence interval, and Wald p-value derived from it,
  was biased by roughly 1.6% in the bundled example. Estimates, log-likelihoods,
  and gradients were unaffected.
* Added `tests/testthat/test-analytic-derivatives.R`, which validates the
  analytic gradient and Hessian of `terradish_algorithm()` against `numDeriv`
  for every measurement and conductance model the package ships, and regression
  tests the constrained-nuisance case above.
* Documented the S3 methods for `terradish_directed` and `terradish_hierarchical`
  objects, which previously had no help pages: `?terradish_directed_methods` and
  `?terradish_hierarchical_methods` now explain what `summary()`, `coef()`,
  `vcov()`, `confint()`, `AIC()`, and `plot()` return and how to read them.
* Added the missing `@return` sections to `NewtonRaphsonControl()` and
  `HagerZhangControl()`.
* Expanded the `melip` dataset documentation with grid dimensions, value ranges,
  units, coordinate system, and a runnable example, and corrected the
  getting-started vignette, which described the data as 34 sites (it is 37) and
  altitude as meters (the layer is pre-rescaled to roughly the unit interval).
* Expanded the package-level help (`?terradish`) into a real entry point:
  main functions grouped by role, plus pointers to every vignette.
* Made the vignette closings consistent: all ten now end with a quick-reference
  workflow block, a summary-of-key-functions table, and a see-also list, and
  `large-landscapes` now uses the same heading levels as the rest.
* Moved the unused generalized-Wishart benchmark harness out of `R/` to
  `inst/benchmarks/generalized-wishart-timing.R`. It was dead code that called
  `set.seed()` without restoring the RNG and depended on the archived
  `RandomFields` package. Dropped the now-unused `nloptr` and `corMLPE` from
  `Suggests`; `terradish` implements MLPE itself and never called `corMLPE`.
* `inst/CITATION` now derives its version and year from the DESCRIPTION instead
  of hard-coding them, and `Authors@R` records both authors' ORCIDs.

terradish 0.0.44
---------
* Prepared the package for CRAN submission: corrected the `LICENSE` file to the
  standard two-field form, xz-compressed the bundled `melip` data, declared the
  `grDevices`, `methods::as`, and `stats::formula`/`stats::optim` imports, and
  extended `.Rbuildignore` to exclude development artifacts.
* Guarded `set.seed()` in `terradish_cv()`, `simulate_covariance_response()`, and
  `simulate.radish()` so they no longer change the user's global random state;
  the RNG is saved and restored on exit.
* `terradish_hierarchical()` now detects when the measurement model profiles to
  `tau = 0` (no detectable resistance-distance signal). The conductance surface
  is unidentified there, so the fit is returned with a warning and `NA` standard
  errors, effective degrees of freedom, and AIC instead of inverting a singular
  Hessian. It also warns when the selected `tau2` sits at the edge of
  `tau2_grid`.
* `terradish_directed()` now fails fast with a clear message when a Wishart
  measurement model is supplied without `nu`.
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
