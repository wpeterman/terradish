# terradish 0.1.0 core release: implementation plan

**Owner:** Bill Peterman. **Prepared:** 24 September 2026. **Status:** ready for handoff. An independent reviewer checked this plan against the source, and its corrections are incorporated.

**Revision (24 September, 11:30):** the owner decided to keep the environmental kernel (`wishart_covariates`) for comparability with MLPE. It stays on `master`, reformulated so that each environmental term is the covariance form of the same pairwise environmental distance that `mlpe_covariates` uses (F8, I11, I12, C10).

**Repository:** `C:\Users\peterman.73\OneDrive - The Ohio State University\R\Packages\terradish` (origin `github.com/wpeterman/terradish`).

**Base commit:** `6a2c13f` ("Release terradish 0.0.47", `master`).

**Evidence:**
- **Audit summary:** a feature audit of 0.0.47, summarized in `review_work_20260924/package_audit/terradish_package_audit_20260924.md`.
- **Reproducing scripts:** `review_work_20260924/package_audit/{core,scale,flex,ext,cv,compute,identifiability}/`, under `C:\Users\peterman.73\OneDrive - The Ohio State University\Research\R_packages\terradish\`.
- **Paths:** the scripts were written in a Linux workspace. Edit paths such as `/tmp/claude-0/...` and `/home/claude/...` before running them.
- **Line numbers** refer to `6a2c13f` and will drift.

**Goal:** release a smaller package whose every advertised feature is sound, tested and documented honestly.
- Development continues on an `experimental` branch that keeps the full 0.0.47 feature set.
- The manuscript will describe only the 0.1.0 scope.

---

## 0. Ground rules for the worker

1. **Follow `AGENTS.md`.**
   - Make small, reviewable diffs.
   - Document with roxygen2 only; never hand-edit `man/`.
   - Add or update tests for every nontrivial change.
   - Write American English, with no em dashes inside sentences.
   - Never fabricate references.
2. **Numerical changes are intentional here, but none may be silent.** Every change to a default, a return value or a numerical result gets a NEWS entry and a test that pins the new contract.
3. **The owner's uncommitted edit.** `vignettes/model-comparison.Rmd` carries an uncommitted edit by the owner, who commits it before handoff (P0.1).
   - If it is still uncommitted when you start, stop and ask.
   - Never stash, reset or overwrite it.
4. **Read-only git commands:** use `git --no-optional-locks ...`, so no stale `.git/index.lock` files are left in the OneDrive folder.
5. **Commit authority.**
   - Make one commit per phase, on the branches named in Section 2, once that phase's acceptance checks pass.
   - Do not push, and do not merge into `master`, without owner sign-off.
   - *(Owner: confirm or edit this item before handoff. Otherwise AGENTS.md forbids commits until asked.)*
6. **Versioning.**
   - Each phase commit bumps the development version by 0.0.1, following the sequence in Section 2. Update DESCRIPTION and NEWS each time.
   - The release commit sets `Version: 0.1.0`.
7. **Ambiguity.** When a statistical choice is not specified here, stop and ask (AGENTS.md). Items marked **Owner decision** give a default; follow it unless the owner changes it before handoff.
8. **Tolerances.** "Unchanged" means:
   - relative differences of at most 1e-6 in logLik;
   - absolute differences of at most 1e-5 in θ, φ and SEs;
   - iteration counts are ignored, and bit-identical results are never required.

---

## 1. Release scope

### 1.1 What 0.1.0 is

Reuse this statement in the README, DESCRIPTION and the package help page.

terradish estimates landscape conductance surfaces from genetic data by maximum likelihood on a raster graph. It fits log-linear conductance models and can optionally estimate each covariate's Gaussian smoothing scale and spline shape. These models are fitted to pairwise genetic distances through the MLPE likelihood, or to allele-frequency covariance (or squared distances derived from it) through a Wishart likelihood on site contrasts. In both families the same pairwise environmental covariates can enter alongside resistance, so isolation-by-environment and isolation-by-resistance estimates are comparable across them. Exact sparse Cholesky and algebraic multigrid solvers handle large rasters. Fixed-domain spatial cross-validation compares conductance formulas.

Estimated coefficients describe relative conductance, and they are conditional on the measurement-model terms. For Wishart likelihoods, the degrees of freedom ν set the precision of every standard error, test and information criterion. Setting ν to the SNP count overstates the information by roughly an order of magnitude.

### 1.2 Feature decisions

- **KEEP:** the feature stays on `master` as it is, apart from documentation.
- **FIX:** the feature stays on `master` after the listed fixes.
- **MOVE:** the feature is removed from `master` and kept on `experimental`.

| Area | Exported objects | Decision | Reason (audit reference) |
|---|---|---|---|
| Graph | `conductance_surface`, `crop_to_focal_buffer`, `terradish_distance`, `terradish_algorithm`, and the deprecated forwarding wrappers `radish_distance`, `radish_algorithm` | FIX | Default `directions` contradicts docs; crop needs caveats; broken example (compute §1, §7) |
| Graph | `terradish_grid`, `radish_grid` | FIX; approximations limited to `"none"` | Landmark and coarse screening inside `terradish_grid` have no exact refinement, so they move to experimental |
| Conductance | `loglinear_conductance` | KEEP | Verified (core §6) |
| Conductance | `gaussian_smoothed_loglinear_conductance`, `gaussian_scale_summary` | FIX | Even-grid misregistration, kernel truncation, NA bounds, σ inference (scale §1–2) |
| Conductance | `smooth_loglinear_conductance` | FIX | Knots are rebuilt in prediction and plotting; reference level of bands (flex §1) |
| Conductance | `linear_conductance` | MOVE | The scale of θ is not identifiable (core §7) |
| Fitting | `terradish`, `radish` | FIX | Convergence reporting, stopping rule, boundary starts, φ warm start (core §6) |
| Scale search | `terradish_scale_optim`, `terra_radish_scale_optim` | FIX | Units, K, default kernel (scale §3) |
| Scale search | `terradish_multiscale`, `radish_multiscale` | MOVE | A warm start that duplicates `approximation = "coarse_raster"` (scale §4) |
| Measurement | `mlpe`, `leastsquares`, `generalized_wishart`, `check_distance_response` | KEEP with small fixes | Verified against brute-force likelihoods (core §2, §4, §5) |
| Measurement | `wishart_covariance` | FIX (reimplement) | Full-rank likelihood on an uncentered E; collapses when sites are clustered (core §1) |
| Measurement | `mlpe_covariates`, `pairwise_endpoint_covariates` | FIX | anova nesting, input validation, README example (ext §3) |
| Measurement | `wishart_covariates` | FIX (reformulate) | **Owner decision (24 Sept): keep, for comparability with MLPE.** The likelihood matched brute force (ext check1). Reformulate each term as the covariance form of the same pairwise environmental distance `mlpe_covariates` uses; contrast form for covariance input; boundary inference for λ; IBE:IBR ratio reporting (F8, I11, I12, C10). |
| Measurement | `wishart_drift_covariates` | MOVE | The diagonal term is confounded with local isolation, sampling variance and rare-variant weighting (ext §2; identifiability addendum) |
| Measurement | `pair_subset_measurement_model` | MOVE | Niche; drops γ when wrapping `mlpe_covariates`; no CV support (ext §4b) |
| Models | `terradish_hierarchical`, `conductance_field` | MOVE | Ignores `conductance_model`; memory; biased conditional AIC; identified only as a generic residual (flex §2; addendum) |
| Models | `terradish_directed`, `terradish_directed_algorithm`, `directed_rates`, re-exports `edge_gradient`, `edge_flow` | MOVE | Not identifiable with gradient covariates; untested otherwise |
| Post-fit | `conductance` method, `mlpe_response_change`, `aic_table`, `slim_terradish`, S3 methods | FIX | Missing `vcov`/`confint`/`nobs`; the `support = "focal"` bug; `x` is ignored (core §8–9; compute §8) |
| CV | `terradish_folds`, `terradish_cv_folds` | FIX; becomes the only CV path | Baseline, failed folds, checkpoints; per-formula conductance models (cv §2) |
| CV (legacy) | `terradish_cv`, `terradish_cv_replicates` (and its methods), `radish_cv`, `cv_model_selection` | MOVE | External-θ scoring bug; φ re-profiling favors measurement extensions; a second CV system (cv §1) |
| Results I/O | `terradish_results`, `terradish_parameters`, `radish_parameters` | MOVE | Legacy results-directory helpers; unaudited |
| Design | `simulate_covariance_response`, `covariance_response_power` | FIX | σ units, `nu_fit`, caveats (cv §5; scale §5) |
| Design | `terradish_assess_settings`, `terradish_solver_benchmark` | KEEP (documentation only) | Speed probes |
| Large graphs | Solvers `direct`, `amg`, `auto` | FIX (the `auto` rule) | compute §2 |
| Large graphs | Solvers `pcg`, `pcg_jacobi`, `block_cg` | MOVE | `block_cg` fails with six or more right-hand sides; the others are experimental |
| Large graphs | `terradish_kron_reduce`, `terradish_kron_reduce_tiled` | MOVE | A diagnostic not used in fitting; falsely described as out-of-core (compute §5) |
| Large graphs | `approximation = "landmark"` and `"coarse_raster"` in `terradish()` | FIX | Keep them only with exact refinement; remove `exact_refine = FALSE` (compute §7) |
| Utilities | `scale_covariates`, `scale_to_0_1`, `lower`, `pca_dist`, optimizer controls | FIX / KEEP | Multi-layer bug in `scale_to_0_1`; `scale_covariates` needs a way to apply stored scales (compute §3) |
| Re-exports | `cov_from_biallelic`, `cov_from_genetic_data`, `fst_from_biallelic`, `dist_from_cov`, `dist_from_biallelic` | KEEP with documentation caveats | landgraph changes in Section 14 |
| Parallel | `cores` argument | KEEP; marked experimental in docs | compute §4 |

### 1.3 New or completed features in 0.1.0

| ID | Feature | Priority |
|---|---|---|
| N1 | `vcov()`, `confint()` and `nobs()` methods for `terradish` fits, including slim fits | P0 |
| N2 | `wishart_covariance()` evaluated on site contrasts, identical to `generalized_wishart(dist_from_cov(C))` | P0 |
| N3 | `terradish_rescale_nu(fit, nu)`, a ν-sensitivity tool for Wishart fits that needs no refit | P0 |
| N4 | `terradish_cv_folds()`: a per-formula `conductance_model` list, and paired fold differences with SEs | P0 |
| N5 | A convergence record in every fit, printed, and flagged in `aic_table()` | P0 |
| N6 | `conductance(new_surface, fit)`: prediction to a new raster for the three core conductance models | P1 |
| N7 | `pairwise_covariates()`: supported input for user pairwise matrices such as geographic distance | P1 |
| N8 | `covariance_response_power(..., nu_fit = )` | P1 |
| N9 | Spline monotonicity report in `summary()` | P1 |
| N10 | Automatic restart when a fit ends on the no-structure boundary | P1 |
| N11 | `gaussian_scale_profile()`: a profile-likelihood interval for σ | P2 |
| N12 | Repeated k-fold support in `terradish_cv_folds()` | P2 |
| N13 | `terradish_cv_folds(..., nuisance = "fixed")`: predictive scoring that holds measurement parameters at training values, so IBR and IBR + IBE models can be compared by CV in either family | P1 |
| N14 | IBE:IBR ratio (γ/β for MLPE, λ/τ for Wishart) with delta-method SE in `summary()`, and an extractor `terradish_ibe_ratio()` | P0 |
| N15 | `wishart_covariates()` built from the same pairwise transforms and objects as `mlpe_covariates()` | P0 |

---

## 2. Branch, commit and sync strategy

The existing `dev` branch (local and `origin/dev`) is stale: it has no commits beyond `master` and is 18 behind. Its name also collides with the `dev/` directory. Leave it untouched.

**Branches.**
1. Run `git tag v0.0.47 6a2c13f`, if the tag does not already exist.
2. After P0.1, create `release/0.1.0` from `master`. All phases happen on this branch.
3. Create `experimental` from `master` at the same point. It keeps every 0.0.47 feature. Add `EXPERIMENTAL.md` (Section 13).

**Commit sequence on `release/0.1.0`:**

| Commit | Content | Version | Portable to `experimental`? |
|---|---|---|---|
| C0 | Phase 0 baseline script and log | 0.0.48 | Yes (normal merge) |
| C1a | Phase 1a refactor only: relocate helpers; no removals and no behavior change | 0.0.49 | Yes (normal merge) |
| C1b | **The split** (Phase 1b): every removal of moved features, in code, tests, roxygen text, vignettes, README and benchmarks, and nothing else | 0.0.50 | **No** (see below) |
| C2–C7 | Phases 2–7: fixes, additions and rewrites of kept features only | 0.0.51–0.0.56 | Yes (normal merge; resolve conflicts in favor of keeping experimental features) |
| Release | Phase 8 | 0.1.0 | Yes |

**Syncing `experimental`** (after the release, or whenever it is wanted):
```
git checkout experimental
git merge <C0> <C1a>                       # normal merges
git merge -s ours <C1b> -m "Record 0.1.0 split without removing experimental features"
git merge master                           # brings C2 onward
```
- **Record C1b first.** Never merge `master` into `experimental` without first recording C1b with `-s ours`; otherwise the deletions propagate.
- **Moved features.** Fixes to moved features happen only on `experimental`.
- **After the merge,** review the package overview and README on `experimental` so they describe its extra features.

---

## 3. Phase 0: preparation and numerical baseline (commit C0)

**P0.1** (owner). Commit the pending `vignettes/model-comparison.Rmd` edit on `master`.
- *Acceptance:* `git --no-optional-locks status` is clean.

**P0.2.** Create the tag and branches from Section 2.
- *Acceptance:* the branches exist.

**P0.3.** Write `dev/release_0.1.0/core_baseline.R` and run it on the base code.
- Put the script in `dev/release_0.1.0/`, which is tracked; `dev/check/` is git-ignored.
- It fits the reference models below with **explicit** settings: `directions = 4`, `curvature = "exact"`, and fixed seeds.
- It saves to `dev/release_0.1.0/core_baseline_0.0.47.rds`: coefficients, logLik, SEs, φ, iterations and the final gradient.
- Models:
  - (a) melip `loglinear_conductance` + `mlpe`, `~ altitude + forestcover`;
  - (b) the same with `leastsquares`;
  - (c) `generalized_wishart` on `dist_from_cov` of a seeded `simulate_covariance_response` draw;
  - (d) `wishart_covariance` on the same draw, uncentered;
  - (e) `smooth_loglinear_conductance` with `s(altitude, df = 4)`;
  - (f) `gaussian_smoothed_loglinear_conductance` on an **odd**-dimension crop of the melip aggregate, with `sigma_upper` at or below one sixth of the smaller raster dimension (see F1);
  - (g) the same on the even-dimension aggregate;
  - (h) `mlpe_covariates` with absolute-difference IBE terms;
  - (i) `wishart_covariates(..., model = "generalized_wishart")` with its current defaults on the draw from (c), using one site covariate.
- *Acceptance:* the RDS file is saved.

**P0.4.** Run `devtools::test()` and `devtools::check()` on the base.
- Record the results, timing and NOTEs in `dev/release_0.1.0/log.md`, together with the final max|gradient| of every fit.
- *Acceptance:* recorded.

---

## 4. Phase 1a: refactor-only relocation (commit C1a)

The goal: no helper that kept code needs may remain in a file that Phase 1b deletes. Behavior must not change, and the test suite must pass unchanged.

**A1.** From `R/radish_multiscale.R`, move **all four internal helpers** that are defined before the `terradish_multiscale` roxygen block into `R/radish_grid.R`:
- `.validate_multiscale_covariates`
- `.aggregate_covariates`
- `.normalize_coarse_raster_control` (holds the `exact_refine` default)
- `.coarse_raster_surface`

Also move `#' @importFrom terra aggregate` (line 163) to `R/terradish-package.R`.

**A2.** Create `R/aic_table.R`. Move `aic_table` into it from `R/radish_cv.R`, with every helper it calls. That includes at least `.default_model_name`, `.model_formula_label`, `.model_rhs_label`, `.mlpe_covariate_count` and the comparison-contract helpers.

**A3.** Create `R/cv_helpers.R`. Move into it the helpers from `R/radish_cv.R` that `R/spatial_cv.R` uses. That includes at least `.terradish_measurement_model` (called at spatial_cv.R:266) and `.fit_terradish_with_fallback` (called at 315, 316 and 321).

**A4.** `.pair_subset_symm` (defined in `R/pair_subset_measurement.R`) is used by the kept `R/wishart_covariates.R` (lines 248, 252, 475 and 481). Move it to `R/wishart_covariates.R`, or to a small `R/utils_matrix.R`.

**A5.** Transitive check.
- After `devtools::load_all()`, run:
  ```
  fns <- ls(asNamespace("terradish"), all.names = TRUE)
  for (f in fns) {
    obj <- get(f, asNamespace("terradish"))
    if (is.function(obj)) {
      g <- codetools::findGlobals(obj, merge = FALSE)$functions
      miss <- g[!vapply(g, exists, logical(1), envir = asNamespace("terradish"))]
      if (length(miss)) cat(f, ":", miss, "\n")
    }
  }
  ```
- It should report nothing new relative to the base.
- Save the script as `dev/release_0.1.0/check_globals.R`; it is rerun after Phase 1b.

**Acceptance.** Test and check results are unchanged from P0.4, and the Phase 0 baseline is unchanged.

---

## 5. Phase 1b: the split (commit C1b; removals only)

For each MOVE item:
1. Delete or edit the listed code, tests, roxygen text, vignette content and README content.
2. If C++ changed, run `Rcpp::compileAttributes()`. Then run `devtools::document()`.
3. Confirm that no Rd files are orphaned (`git status`).

Vignette and README content that *uses* a moved feature is removed in this phase, so every check from here on can build the vignettes. Phase 7 rewrites the remaining kept content. Where a removal leaves a gap, insert an HTML comment (`<!-- Phase 7: rewrite -->`), never placeholder prose.

### S1. Directed model
- **Code:**
  - Delete `R/directed_conductance.R` and `src/directed_sparse_lu.cpp`.
  - Remove the `edge_gradient` and `edge_flow` blocks from `R/reexports.R`.
  - Regenerate the Rcpp exports.
- **Tests:** delete `test-directed-conductance.R`.
- **Docs, vignettes and README:**
  - Delete `vignettes/directional-conductance.Rmd` and `vignettes/vignette-directional.rds`.
  - Remove the `vignette("directional-conductance")` links: getting-started (around lines 887 and 890), large-landscapes (around 140 and 459), and any others grep finds.
  - Remove the directed mention in large-landscapes.
  - Edit `R/terradish-package.R` at line 41 and lines 52–53.
  - Fix the header comment at `vignettes/precompute.R:9`.
  - Delete `dev/make_vignette_directional.R`.

### S2. Hierarchical field
- **Code:** delete `R/hierarchical_conductance.R`.
- **Tests:**
  - Delete `test-hierarchical-conductance.R`.
  - In `test-verbose-and-tau2.R`, remove the two hierarchical blocks (lines 91 and 149) and keep the three optimizer tests.
- **Docs, vignettes and README:**
  - Delete `hierarchical-conductance.Rmd` and `vignette-hierarchical.rds`.
  - Remove the `vignette("hierarchical-conductance")` links in getting-started, spline-conductance (around line 983), wishart-covariance (around 845), ibe-ibr (around 732) and large-landscapes.
  - Edit `R/terradish-package.R` at lines 39–40 and 52–53.
  - Delete `dev/make_vignette_hier.R`.

### S3. Drift covariates
- **Code:** delete `R/wishart_drift_covariates.R`.
- **Tests:**
  - Delete `test-wishart-drift-covariates.R`.
  - In `test-analytic-derivatives.R`, remove only the drift parts of the mixed blocks (around lines 116, 168 and 195). Keep their spline and `mlpe_covariates` parts.
- **Docs, vignettes and README:**
  - Remove the `wishart_drift_covariates` entries from the `@seealso` lines at `R/generalized_wishart.R:106-107`, `R/wishart_covariance.R:120-121` and `R/wishart_covariates.R:89` (keep `wishart_covariates`).
  - In `wishart-covariance.Rmd`, delete "Modeling site-specific diagonal variance" (lines 523–637, including the `drift-fit` chunk).

### S4. Environmental kernel: retained
Owner decision (24 September): `wishart_covariates` stays on `master` and is reformulated in Phases 2–4 (F8, I11, I12, C10). Phase 1b removes only its pair-subset hooks (S6) and its uses of the legacy CV (S10). The kernel blocks of `test-analytic-derivatives.R` (around line 78) stay; they also keep the coverage for a profile Hessian with a nuisance parameter at its bound.

### S5. `linear_conductance`
- **Code:**
  - Delete its roxygen block and function in `R/radish_conductance_model.R`, from the documentation preceding line 597 through the class assignment at line 652.
  - Edit the mentions at lines 272, 328 and 343.
  - Edit the class allow-lists at `R/radish_grid.R:284` and `:701`.
- **Tests:** remove the linear parts of `test-analytic-derivatives.R`.
- **Docs:** edit `R/terradish-package.R:28` and the getting-started mention (around line 855).

### S6. Pair subsets
- **Code:**
  - Delete `R/pair_subset_measurement.R`.
  - Remove the hooks at `R/radish_optimize.R:258`, 265, 325–331 and 1430.
  - Remove the hooks in `R/aic_table.R` (formerly radish_cv.R:159–182 and 276–283). The BIC pair count reverts to n(n−1)/2.
  - **Keep** the `subsetter` in `R/mlpe_covariates.R:411`; CV uses it.
- **Tests:**
  - Delete `test-pair-subset-measurement.R`.
  - Remove the pair-subset blocks from `test-model-selection.R` and `test-analytic-derivatives.R`.
- **Docs:** remove the README section and the getting-started mention.

### S7. Multiscale wrapper
- **Code:** delete `R/radish_multiscale.R`. Its helpers were relocated in A1, and `radish_multiscale` is defined in this same file.
- **Tests:** remove `radish_multiscale` from `test-legacy-wrappers.R`.
- **Docs:** remove any mention in the vignettes or README.

### S8. Kron reduction
- **Code:** delete `R/kron_reduction.R`.
- **Tests:** delete `test-kron-reduction.R`.
- **Docs:** in `large-landscapes.Rmd`:
  - delete "An exact, memory-bounded reduction" (lines 344–410, including the `reduce` chunk);
  - delete the Kron line in the quick reference (around line 432);
  - delete every "out-of-core" phrase (lines 67–68 and 404).
  - NEWS history stays as it is.

### S9. Experimental solvers
- **Code.**
  - In `R/radish_algorithm.R`, remove `pcg`, `pcg_jacobi` and `block_cg` from:
    - the `match.arg` at line 333;
    - the defaults at lines 303–305;
    - the dispatch at lines 568–570;
    - the `@param solver` text at line 707;
    - the exported `terradish_algorithm()` signature at line 765.
  - Remove them from `R/radish_optimize.R:917`, from the docs at `R/radish_optimize.R:690`, and from `R/covariance_response_power.R:168`.
  - **C++:** if the solver functions share `src/radish.cpp` with kept code, leave them compiled but unexported. Otherwise remove them and regenerate the exports.
- **Tests:** delete `test-block-cg.R`.
- **Docs:**
  - Remove the `block_cg` and `pcg_jacobi` mentions in `large-landscapes.Rmd`.
  - Update `inst/benchmarks/synthetic-solver-scaling.R` to use `direct` and `amg` only.

### S10. Legacy CV and results I/O
- **Code:**
  - From `R/radish_cv.R`, remove:
    - `terradish_cv`;
    - `terradish_cv_replicates` and its print and summary methods;
    - `radish_cv` and `cv_model_selection`;
    - `terradish_results`, `terradish_parameters` and `radish_parameters`;
    - the private helpers used **only** by those functions.
  - Everything kept was already relocated in A2 and A3. If the file ends up empty, delete it.
- **Tests:**
  - Delete `test-cv-helpers.R`, which is entirely legacy.
  - Remove the `cv_model_selection` blocks in `test-model-selection.R` (lines 37–62 and 94–100) and in `test-spatial-cv.R`.
  - Remove the `terradish_cv(` block in `test-ibe-workflow.R` (around line 112).
  - Remove the `terradish_cv` blocks in `test-landmark-grid-cv.R` (lines 67 and 118); `terradish_cv_folds` has no equivalent.
  - Remove `radish_parameters` and `radish_cv` from `test-legacy-wrappers.R`.
- **Docs, vignettes and README:**
  - Replace the `\link{terradish_cv}` links with `\link{terradish_cv_folds}` in `R/mlpe_covariates.R` (lines 50, 62, 147, 152 and 155) and `R/terradish-package.R:44`.
  - In `model-comparison.Rmd`, delete the legacy CV sections (lines 481–632) and leave a Phase 7 comment. Preserve the owner's P0.1 text elsewhere in the file.
  - In `ibe-ibr-workflow.Rmd`, delete the CV section (lines 313–379) and "Step 7: Saved-results utilities" (lines 437–476).
  - Delete the README's legacy CV sections.

### S11. Unrefined approximations
- **Code:**
  - Make `exact_refine = FALSE` an error in `terradish()` and its helpers: "not supported in terradish 0.1.0; available on the experimental branch".
  - `terradish_grid()` accepts only `approximation = "none"`.
- **Tests:** remove the tests that use `exact_refine = FALSE` or approximations in `terradish_grid`.
- **Docs:** remove the `exact_refine = FALSE` text at README:287 and in the `approximation_control` docs (`R/radish_optimize.R:716-743`).

### Acceptance for Phase 1b

1. `devtools::document()` runs clean, and rerunning `dev/release_0.1.0/check_globals.R` reports nothing.
2. This grep over `R/ src/ tests/ vignettes/ inst/ README.md` returns nothing:
   ```
   grep -rnE "terradish_directed|directed_rates|edge_gradient|edge_flow|terradish_hierarchical|conductance_field|wishart_drift_covariates|(^|[^g_])linear_conductance|pair_subset|terradish_multiscale|radish_multiscale|terradish_kron_reduce|kron_reduce|block_cg|pcg_jacobi|\"pcg\"|terradish_cv[^_]|terradish_cv_replicates|radish_cv|cv_model_selection|terradish_results|terradish_parameters|radish_parameters|hierarchical-conductance|directional-conductance|out-of-core|exact_refine = FALSE"
   ```
3. `devtools::test()` passes.
4. `devtools::check()`, **with vignettes built**, shows no new errors or warnings relative to P0.4.
5. Baseline models (a), (b), (c), (e), (f) and (h) are unchanged, within the Section 0 tolerances.

---

## 6. Phase 2: correctness fixes (commit C2)

### F1. Gaussian crop and kernel support
**Location:** `R/gaussian_scale_conductance.R:224-227` (`.gaussian_scale_prepare_layer`), and every use of `i1:i2` or `j1:j2`. Also grep for `floor(nr / 2)` and `floor(nc / 2)`.

**(a) Crop.**
- Set `i1 = ceiling(nr/2)` and `i2 = ceiling(nr/2) + nr - 1`, and do the same for columns. Odd dimensions stay unchanged.
- The reviewer confirmed that with this fix the smoother matches a brute-force convolution over the same support to 1e-15. Without it, 8 × 10 grids were off by 0.18–1.9.

**(b) Kernel support.**
- The problem: the kernel matrix is nr × nc, so offsets beyond about ±n/2 are truncated, asymmetrically on even grids.
- **Owner decision; default:** keep the truncated kernel, but set the default `sigma_upper` to (smaller raster dimension in cells)/6 × cell size, so the window spans at least ±3σ.
- Warn when a user sets `sigma_upper` above that value, and document this in the help page.
- The alternative, a full symmetric (2n−1) kernel, would change results at large σ on every grid and needs a separate decision.

**Tests** (`test-gaussian-alignment.R`):
1. On 7 × 9 and 8 × 10 rasters with an NA mask, for σ ∈ {0.7, 1.5, 3} cells, the package smooth equals a brute-force normalized convolution **over the same kernel window** to 1e-10.
2. A spike at (4, 5) on an even raster peaks at (4, 5).
3. The σ derivatives match numDeriv on even grids.
4. The default `sigma_upper` rule applies, and the warning fires above it.

### F2. `wishart_covariance` on contrasts (N2)
**Location:** `R/wishart_covariance.R`, the objective (lines 185–215) and its derivatives.

**Specification.**
- Let L be an orthonormal n × (n−1) contrast basis. The objective is (ν/2)[log det Σ̃ + tr(Σ̃⁻¹S̃)], with Σ̃ = τLᵀEL + σI and S̃ = LᵀSL.
- Map derivatives back with ∂/∂E = L(∂/∂Ẽ)Lᵀ.
- Keep the φ names and their meanings.
- **Do not** add a ν ≥ n−1 rule; `generalized_wishart` has none.
- `fitted` stays Σ = τE + σI (n × n). Residuals are H(S − Σ)H.
- Mirror the structure of `R/generalized_wishart.R`.
- Before coding, establish how σ maps between the two parameterizations. GW's nugget may enter the distance model as 2σ.

**Tests** (rewrite `test-wishart-covariance.R`):
1. For a site-centered C, logLik, θ̂ and the fitted contrasts equal those from `generalized_wishart(dist_from_cov(C))` to 1e-6.
2. Adding a·11ᵀ to C, or centering C, leaves the fit unchanged.
3. Regression: 20 sites in one corner of a 30 × 30 raster (audit script `core/exp1b_wishart_clustered.R`), three seeds. Require τ̂ > 0 and θ̂ within 3 SE of the truth.
4. The analytic gradient and Hessian match numDeriv.

### F3. Spline knots reused everywhere
**Location:**
- `smooth_loglinear_conductance` (`R/radish_conductance_model.R`, spec handling at lines 100–115);
- `conductance()` with `support` (`R/radish_graph.R:390-401`);
- the plot factory path (`R/plot_terradish.R`);
- the coarse stage (`R/radish_optimize.R`, around line 1191).

**Specification.**
- At fit time, store the full basis specification of each `s()` term in the fit: basis, df, degree, interior and boundary knots, and the centering from F4.
- Every rebuild uses that stored specification: focal support, plotting with a factory, the coarse stage, and prediction (R7).
- Knots are never recomputed from new data.

**Tests** (`test-spline-support.R`):
1. `support = "focal"` equals `support = "none"` on unclamped cells, to 1e-10.
2. The plot curve drawn with a factory equals the curve drawn without one.
3. A coarse warm start followed by exact refinement equals a direct fit.

### F4. Centered spline basis
**Location:** `smooth_loglinear_conductance`.

**Specification.**
- Center each basis column over the active graph cells, and store the centers.
- logLik and θ̂ do not change, because a constant shift in log-conductance is absorbed.
- Bands then refer to the landscape-mean log conductance.
- `conductance()` values shift by a constant; note this in NEWS.

**Tests:** logLik and θ̂ equal baseline (e) within tolerance, and the band width at the covariate minimum is greater than 0.

### F5. Outer-search scale units
**Location:** `R/terradish_scale_optim.R`: the docs at line 229, the example, and the default `scale_fun`.

**Specification.**
1. With the multiScaleR kernel, σ is in **map units**, as the audit verified. Say so in the docs.
2. Make `scale_fun = "terradish"` the new default. It uses the package's own normalized Gaussian (F1), so outer-search estimates match the joint model.
3. Keep `"multiScaleR"` as an option, and document that its σ is not interchangeable with terradish's.
4. Count the scale parameters in K.
5. Fix the bounds in the melip example.

**Tests:**
1. The default path reproduces the joint σ̂ to within 1%.
2. Its df equals the joint model's df.
3. The multiScaleR path uses map units on a resolution-100 raster.

### F6. σ units in the simulator
**Location:** `simulate_covariance_response` (`R/genetic_distance.R`).

**Specification:** accept `sigma.<layer>` in the map units that `coef()` reports, and convert internally.

**Test:** a round trip on a resolution-100 raster.

### F7. `scale_to_0_1` per layer
**Location:** `R/convenience_funcs.R:45-53`.

**Specification:** use a separate range for each layer, as documented.

**Test:** a two-layer raster.

### F8. Environmental kernel reformulation (N15; owner decision to keep)
**Location:** `R/wishart_covariates.R`: the constructor (lines 134–164), `.make_wishart_kernel_covariates` (166–205), the covariance-path objective (near line 332) and the boundary flag (491).

**Why.** The likelihood is correct: it matched a brute-force dense computation for both paths (audit `ext/check1`). Three things limit its validity and its comparability with MLPE:
1. The covariance path uses the full-rank likelihood, with the same centering flaw as `wishart_covariance` (F2).
2. Each kernel is the outer product z_c z_cᵀ. In distance form that is λ(z_i − z_j)², the `sqdiff` transform, whereas `mlpe_covariates` defaults to `absdiff`.
3. `normalize = TRUE` divides by mean(z_c²). That makes `scale` inert and puts λ on a different scale from the MLPE coefficient.

**Specification.**
- **(a) Contrast form.** Evaluate the covariance path on site contrasts exactly as in F2. It must equal the `model = "generalized_wishart"` path given `dist_from_cov(S)`.
- **(b) Kernels from the MLPE pairwise transforms.** The new signature mirrors `pairwise_endpoint_covariates()`:
  `wishart_covariates(x, coords = NULL, transform = c("absdiff", "sqdiff", "euclidean", "manhattan"), scale = FALSE, model = c("wishart_covariance", "generalized_wishart"))`.
  - `x` may also be a `terradish_pairwise_covariates` object, from `pairwise_endpoint_covariates()` or `pairwise_covariates()` (R8). One object then feeds both `mlpe_covariates()` and `wishart_covariates()`.
  - For each pairwise covariate k, build the n × n distance D_k and the kernel K_k = −½ H D_k H, where H is the centering matrix.
  - **PSD check.** All four transforms are conditionally negative definite, so each K_k is positive semidefinite; I checked this numerically on random inputs (smallest eigenvalue ≥ −1e−15 × largest). Still verify every kernel (smallest eigenvalue ≥ −1e−10 × largest) and error otherwise; user matrices supplied through `pairwise_covariates()` can fail.
  - **Shared additive structure.** Because K_ii + K_jj − 2K_ij = D_ij exactly, the Wishart model's expected squared distance is τR_ij + Σ_k λ_k D_k,ij + 2σ. That is the same additive structure as the MLPE mean α + βR_ij + Σ_k γ_k D_k,ij.
  - `transform = "sqdiff"` reproduces the current kernel z_c z_cᵀ exactly.
  - **Default `transform = "absdiff"`**, matching `pairwise_endpoint_covariates()`. This changes the default kernel (breaking; NEWS).
- **(c) Scaling.**
  - Remove `normalize`.
  - `scale = TRUE` standardizes site values before the transform, exactly as in `pairwise_endpoint_covariates()`. λ_k is then per unit of the same pairwise environmental difference that the MLPE γ_k uses.
  - Recommend `scale = TRUE` in docs and examples.
- **(d) Constraint.** A valid covariance requires λ_k ≥ 0, whereas MLPE γ_k is unconstrained. Document that a negative association (greater genetic similarity with greater environmental difference) appears as λ_k = 0 in the Wishart family.
- **(e) Subsetter.** The CV subsetter slices D_k to the retained sites and rebuilds K_k from that block. It never recomputes transforms from site values standardized within the subset.

**Tests** (rewrite `test-wishart-covariates.R`):
1. Both paths match a brute-force dense likelihood (adapt `ext/check1`).
2. The covariance path equals the generalized path on `dist_from_cov(S)`, to 1e-6.
3. `sqdiff` reproduces z_c z_cᵀ; every transform gives a PSD kernel; a user matrix that is not conditionally negative definite errors.
4. λ recovery: simulate S ~ W(ν, τE + λK + σI)/ν at the true ν on a small graph, with three seeds. λ̂ falls within 3 SE of the truth.
5. One `pairwise_endpoint_covariates()` object fits in both `mlpe_covariates()` and `wishart_covariates()`. R8 later adds the same test for `pairwise_covariates()` objects.
6. Subsetter: a fit on a site subset equals a fit built directly from that subset's data.
7. Baseline (i) refit with `transform = "sqdiff"`: logLik and θ are unchanged, and λ_new = λ_old / mean(z_c²).

### Acceptance for Phase 2
- All tests pass.
- Against the baseline:
  - (d) and (g) change by design;
  - (e) keeps θ and logLik;
  - (i) matches the F8 test 7 relation;
  - all other models are unchanged.
- Record the before and after values in `log.md`.

---

## 7. Phase 3: inference, reporting and ν (commit C3)

### I1. `vcov`, `confint` and `nobs` (N1)

**Specification:**
- **`vcov`:** the inverse of the negative Hessian (θ and any scale parameters), taken from the stored curvature.
- **`confint`:** Wald intervals for θ. For Gaussian σ, give a Wald interval truncated at the bounds, with a note. Use the N11 profile interval when it is available.
- **`nobs` convention:** return the number of pairs that `aic_table()` uses for BIC (n(n−1)/2). Then `stats::BIC(fit)` matches `aic_table(BIC = TRUE)`.
  - AICc keeps its site-count n (radish_cv.R:253–257).
  - Document both conventions, and explain why neither is a true independent sample size.
- Keep the Hessian in slim fits, and register the S3 methods.

**Tests:**
- Values match `summary()`.
- Slim fits give the same results as full fits.
- `BIC(fit)` equals the table value.

### I2. Convergence record (N5)

**Specification:**
- Where the fit object is assembled (`R/radish_optimize.R`, around lines 1345–1377), store `fit$convergence = list(code, message, iterations, max_abs_projected_gradient, criterion, boundary, restarted)`.
- Show one line in `print()` and `summary()`.
- `aic_table()` gains `converged` and `boundary` columns, and warns when any row fails either check.

**Test:** a fit with `maxit = 2` is flagged.

### I3. Stopping rule

**Current behavior:** `R/newton_raphson.R:186` and `R/bfgs.R:82` stop when **either** the gradient test or ftol is met. The docs (`newton_raphson.R:30-31`) say both must be met.

**New rule:**
- Use the **projected** gradient: `gradient_box`, which zeroes the components at active bounds. It is currently computed at line 190 (newton_raphson.R) and line 86 (bfgs.R); compute it before the convergence check.
- Code 0 when max|projected gradient| < ctol.
- Code 0 also when ftol is met **and** max|projected gradient| < sqrt(ctol).
- Code 2 ("stalled") when ftol is met but the projected gradient is larger than that; also issue a warning.
- Update the docs to match.

**Validate against P0.4.** No baseline fit may become code 2. That includes the flagship melip MLPE fit, which stopped at max|gradient| = 1.61e-8. If any does, stop and report to the owner.

**Tests:** a toy objective that stalls returns code 2, and every Phase 0 baseline fit returns code 0.

### I4. Boundary restart (N10)

**Specification:**
- Trigger: the final fit lies on the no-structure boundary (MLPE β = 0 or Wishart τ = 0; see the gradient zeroing at `R/radish_algorithm.R:966`), and the fit did not start from the default.
- Refit once from the default start. Reset only the log-linear coefficients to 0; keep any Gaussian σ at its default start.
- Keep whichever fit has the higher logLik, and set `restarted = TRUE`.

**Test:** audit case `core/exp7_optim.R`, where a start of (2, 2) now reaches logLik 269.7 instead of 207.

### I5. φ warm start

**Current behavior:** at `R/radish_subproblem.R:104,116-122`, `phi_start$phi` fails on an atomic vector, and the code falls back to the default start without saying so.

**Specification:**
- Pass a proper list with `phi`, `lower` and `upper`.
- Keep the fallback, but `message()` when it triggers in verbose mode.
- This changes optimization paths, so baseline agreement is judged within tolerance only.

**Test:** the warm start is used, and there is no silent fallback.

### I6. MLPE ρ at the boundary

**Specification:** in `R/mlpe.R`, when logit(ρ) < −8:
- report ρ at its 0 bound;
- set its SE to NA and exclude it from Wald output;
- note this in `summary()`;
- keep df unchanged, and document why.

**Test:** audit script `core/exp9_rho0.R`.

### I7. ν rescaling (N3)

**Specification:** add `terradish_rescale_nu(fit, nu)` for the Wishart families only. The objective is exactly ν·f, so the function:
- multiplies every stored quantity that scales with ν by ν_new/ν_old: `loglik`, `mle$hessian`, `hessian_internal`, `fit$hessian`, `fit$phi_hessian`, and any cached curvature;
- recomputes `aic` and the SEs;
- sets `comparison$nu`, the field `aic_table` checks (`R/radish_optimize.R:334-341`);
- records `nu_original`;
- leaves the estimates unchanged.

**Tests:**
- The result equals a refit at the new ν, within tolerance.
- `aic_table` refuses a table mixing rescaled and original-ν fits.

### I8. anova nesting

**Location:** `R/radish_optimize.R:1733-1768`.

**Specification:**
- Require identical response, family, ν, graph and conductance factory class.
- Require nested measurement-model covariate sets: the smaller model's `mlpe_covariates` columns must be a subset of the larger model's.
- Otherwise raise an error that names the mismatch.
- Warn when a tested parameter's null value lies on a bound.

**Tests:** non-nested sets raise an error (the audit found χ² = −13.5 in this case), and a nested pair works.

### I9. Gaussian σ reporting

**Specification:**
- Fix the `sigma_lower` and `sigma_upper` indexing at `R/gaussian_scale_conductance.R:1262-1263`.
- `summary()` reports σ̂ and its SE with **no z or p**, and flags σ̂ within 1% of a bound.
- `plot(type = "sigma")` truncates at the bound.

**Tests:** the bounds are filled, the flag works, and σ has no p-value.

### I10. Profile interval for σ (N11, P2)

**Specification:**
- Equal bounds currently raise an error (lines 157–159). Add an internal fixed-σ path, for example an internal `sigma_fixed` argument to the factory.
- Then add `gaussian_scale_profile(fit, layer, n = 25, level = 0.95)`.

**Test:** the interval covers the truth in a seeded case.

### I11. Boundary inference for λ
- **Boundary flag.** Include λ_k = 0 in the flag. Only τ is checked today (`R/wishart_covariates.R:491`).
- **anova.** When the tested parameters are λ's with null value 0, use the chi-bar-square reference:
  - one λ: p = ½P(χ²₁ > LR), from the ½χ²₀ + ½χ²₁ mixture;
  - k λ's: mixture weights C(k, j)2⁻ᵏ on χ²_j. Document this as exact when the λ's are information-orthogonal and approximate otherwise.
  - Print which reference was used.
  - MLPE γ_k is unconstrained, so its tests stay χ².
- **summary.** At λ̂_k = 0:
  - flag the estimate as "at bound";
  - give no symmetric Wald interval;
  - report a one-sided upper limit: λ̂ + z₁₋α SE, or the profile limit if implemented.
- **Tests:**
  - p-values match a hand computation for k = 1 and k = 2;
  - the boundary flag is set;
  - no two-sided interval is given at the bound.

### I12. IBE:IBR ratio (N14)
- **What to report.** `summary()` reports γ_k/β for `mlpe_covariates` fits and λ_k/τ for `wishart_covariates` fits.
  - Label it "resistance-distance equivalent of one unit of environmental difference".
  - In both families it is the ratio of the environmental coefficient to the resistance coefficient in the same additive model for pairwise differentiation (F8b). So, given the same transform and scaling, the ratio means the same thing in both.
- **SE.** Use the delta method on the joint covariance of the two parameters. Return NA, with a note, when β or τ is at 0.
- **Extractor.** Export `terradish_ibe_ratio(fit)`, returning a data frame with covariate, ratio, SE, lower, upper and family.
- **Document:**
  - absolute coefficients are not comparable across families, because the responses differ; the ratio is comparable. It is exactly comparable when the MLPE response is the same covariance-derived distance, and approximately so for distance measures that scale linearly with it. Each family estimates θ separately, so R differs slightly between the fits;
  - the same ratio idea underlies BEDASSLE (Bradburd, Ralph and Coop 2013, *Evolution*). Verify the citation in Zotero before citing it.
- **Tests:**
  - the delta-method SE matches a numerical-derivative SE;
  - fitting both families to the same covariance-derived distance, simulated with λ > 0, gives finite ratios of the same sign.

---

## 8. Phase 4: cross-validation (commit C4)

### C1. Per-formula conductance models (N4)
- **Specification:** `terradish_cv_folds(..., conductance_model = )` accepts either one factory or a named list matching `names(formulas)`.
- **Test:** compare `~ x` log-linear, `~ x` Gaussian and `~ s(x)` spline in one call.

### C2. Uniform baseline
- **Specification:**
  - Score the baseline with `loglinear_conductance` at θ = 0, on the same graph and response and under the same measurement model.
  - Build its formula from the first raster variable returned by `all.vars()` on the first formula's right-hand side, because `loglinear_conductance` cannot parse `s()`.
  - Conductance is then uniform, whatever the formulas.
- **Test:** the baseline is finite when every formula is Gaussian (the audit's `cv/e5.R` case).

### C3. Failed folds
- **Location:** `R/spatial_cv.R:394-401`.
- **Specification:**
  - A model with a failed fold gets an `NA` total and a count, `n_folds_ok`.
  - Add `common_folds_total`, summed over the folds in which every model succeeded, and rank on it.
  - Warn, listing the failures.
- **Test:** force a failure.

### C4. Checkpoint safety
- **Location:** `R/spatial_cv.R:277-284`.
- **Specification:**
  - Store a signature with no new dependencies: `tools::md5sum()` of a temporary file holding `serialize()` of a list.
  - That list contains the response matrix, the graph dimensions, the focal cells, the covariate names with per-layer sums, the formula strings, the conductance and measurement model identities, ν and the folds.
  - A signature mismatch on resume is an error.
- **Test:** changing the response raises an error.

### C5. Paired comparison output
- **Specification:** there is no `summary` method today. Extend the `summary` data frame that `terradish_cv_folds` returns, and `print.terradish_cv_folds`:
  - for each model, the per-fold difference from the best model, the SE of the mean difference, and the number of folds;
  - rank by total held-out logLik;
  - flag models within 1 SE of the best.
- **Tests:** a structure test, plus a seeded case with a known ordering.

### C6. One measurement model per call
- **Specification:** when `nuisance = "reprofile"` (the default), assert that each call uses a single measurement model and a single ν (already true by design). Document that re-profiled scores cannot be compared across calls that use different measurement models or ν. C10 adds the exception.
- **Test:** existing tests.

### C7. Repeated k-fold (N12, P2)
- **Specification:** `folds` may be a list of fold vectors. Summaries pool across the repeats and report the between-repeat SD.
- **Test:** two repeats.

### C8. ν invariance regression
- **Specification:** for log-linear models, fold totals at ν = 500 equal 20 times those at ν = 25 (within 1e-6 relative), and the rankings are identical.
- **Test:** passes.


### C10. Predictive CV for measurement-model comparisons (N13)
- **Specification.**
  - `terradish_cv_folds(..., nuisance = c("reprofile", "fixed"))`, and `model` may be a named list matching `names(formulas)`.
  - With `"fixed"`, every φ, including λ_k and the MLPE γ_k, is held at its training estimate. The test block is scored as an out-of-sample predictive density.
  - Only in this mode may the formulas in one call use different measurement models, for example IBR against IBR + IBE. They must still share one likelihood family, response and ν.
- **Tests.**
  - The audit's pure-noise environmental covariate (`cv/e2.R`), 8 seeded replicates, for both `mlpe_covariates` and `wishart_covariates`. Under `"fixed"`, the extension is not systematically preferred (at most 3 of 8); under re-profiling it was preferred in 8 of 8.
  - A simulated real environmental effect (λ > 0) is preferred in most replicates.

### Interpretation, for the Rd in Phase 6
- **What the score is:** the marginal log-likelihood of the test block, with nuisance parameters re-profiled on each test fold.
- **Valid use:** with `nuisance = "reprofile"`, comparing conductance formulas under one measurement model. With `nuisance = "fixed"`, the score is a predictive density with every measurement parameter held at its training value; it can compare measurement models within one family, such as IBR against IBR + IBE.
- **ν:** rankings do not depend on ν, but gains scale with it.
- **Spatial folds:** they score within-cluster pairs only.
- **Train–test pairs:** pairs between training and test sites are never scored.

---

## 9. Phase 5: robustness and usability (commit C5)

**R1. Default `directions = 8`.**
- Change the default in `R/radish_graph.R:177`, matching the recommendation documented at lines 75 and 122. Record it in NEWS as a changed default.
- Tests pinned to the old default pass `directions = 4` explicitly.
- *Test:* the default.

**R2. Graph messages.**
- Warn when two focal sites share a cell.
- Replace the "located on a missing cell" message for sites on pruned components.
- Document the edge weight c_i + c_j and the absence of diagonal distance weighting.
- *Tests:* the messages.

**R3. `auto` solver** (`R/radish_algorithm.R:176-184`).
- Above the large-graph threshold, choose AMG whatever the number of right-hand sides, unless the user sets `solver`.
- Document that AMG iteration counts rise with conductance contrast.
- *Test:* the selection rule.

**R4. `cores`.**
- Document it as experimental:
  - it uses PSOCK on every OS;
  - only Hessian and partial solves run in parallel;
  - it is disabled for AMG and cached CHOLMOD;
  - it is slower on small graphs.
- Issue a one-time `message()` when `cores > 1` and the graph has fewer than 50,000 cells, or when the solver ignores it.
- *Test:* the message.

**R5. Cropping.**
- Fix the `crop_to_focal_buffer` example; `surface$dim` does not exist.
- Issue a `message()` when the buffer is less than half the maximum distance between sites.
- Document the measured bias (resistance +15–23% at a 2-cell buffer, +1–2% at 10 cells) and recommend a sensitivity refit.
- Repeat the caveat in `?conductance_surface` and `?terradish`.
- *Test:* the message.

**R6. Applying stored scaling.**
- `scale_covariates()` already stores its centers and scales in the `"terradish_scale"` attribute (`R/convenience_funcs.R:167`), and `center` and `scale` are already logical arguments.
- Add a `reference =` argument that takes a previous `scale_covariates()` result, or its attribute, and applies those centers and scales to new rasters.
- *Test:* a round trip.

**R7. Prediction to new rasters (N6).**
- `conductance(x, fit)`, with `x` a new `terradish_graph`, uses the stored specification:
  - the log-linear terms;
  - the stored spline specification (F3);
  - the Gaussian σ with the stored post-smoothing standardization.
- It never re-standardizes on the new raster, and it errors for other factories.
- Document that the user scales the new covariates with R6.
- *Tests:* on the fitting graph it reproduces the fitted surface; on a shifted raster it equals a manual computation.

**R8. Pairwise covariates (N7).**
- A new `pairwise_covariates(...)` combines `pairwise_endpoint_covariates` outputs with named n × n matrices or `dist` objects into one classed object. Both `mlpe_covariates()` and `wishart_covariates()` accept it (F8).
- `mlpe_covariates()` raises an error on an unclassed matrix whose row count is not n(n−1)/2.
- Fix the README example: `pairwise_endpoint_covariates(melip.coords, covariates)` has its arguments reversed.
- *Tests:* geographic distance plus a climate difference fit together, in both `mlpe_covariates()` and `wishart_covariates()`; malformed input raises an error.

**R9. `nu_fit` (N8).**
- `covariance_response_power(..., nu_fit = nu)` simulates at `nu`, fits at `nu_fit`, and reports coverage and power.
- *Test:* coverage falls when `nu_fit` is 20 times `nu`.

**R10. Spline monotonicity (N9).**
- For each `s()` term, `summary()` reports whether the fitted curve is monotone over the covariate range at the focal sites, and how many times its derivative changes sign.
- *Tests:* a monotone case and a hump-shaped case.

**R11. Landmark approximation in `terradish()`.**
- Keep it, with exact refinement mandatory (S11).
- *Test:* landmark plus refinement equals a direct fit within tolerance.

**R12. Diagonal mode check.**
- If the input carries landgraph's `diagonal` attribute set to `"within"` (Section 14), warn that the diagonal is on another scale and recommend `"gower"`.
- *Test:* the warning (skip when landgraph is older).

---

## 10. Phase 6: function documentation (commit C6)

Follow the house style (`r-doc-style`):
- plain-language descriptions of inputs, and outputs a reader can interpret;
- a "how to read the output" walkthrough, `\seealso`, and lightweight runnable examples (`\donttest{}` where slow);
- no em dashes, and American English.

### Package-level and fitting topics

- **`?terradish-package`**
  - Use the statement from Section 1.1 and list the kept functions by family.
  - Add a short paragraph on the `experimental` branch, with `remotes::install_github("wpeterman/terradish@experimental")`.
- **`?terradish`**
  - The convergence record, the restart and the stopping rule (I2–I4).
  - The solvers are `direct`, `auto` and `amg`.
  - Approximations require exact refinement.
  - The standard ν paragraph (Appendix C).
  - The crop caveat under "Large rasters".
- **`?terradish_grid`:** only `approximation = "none"`.
- **`?conductance_surface`:** the default `directions = 8`, the edge weight, the crop caveat, and the duplicate-cell warning.

### Conductance models and scale search

- **`?loglinear_conductance`:** coefficients are relative log conductance per SD of each covariate.
- **`?gaussian_smoothed_loglinear_conductance`, `?gaussian_scale_summary`**
  - σ is raster smoothing in map units.
  - β is per SD of the smoothed layer, which is re-standardized at every σ.
  - The half-cell lower bound still smooths.
  - Near the upper bound the smoothed layer acts as a broad trend covariate.
  - Dispersal averaging produces σ > 0 even without a real scale of effect.
  - Document the truncated kernel and the new default `sigma_upper` (F1), the absence of z and p for σ, the bound flag, and the profile interval (I10).
- **`?smooth_loglinear_conductance`**
  - It is an unpenalized regression spline with fixed df; drop "GAM-like" and "partition of unity".
  - The basis is centered, so bands are relative to the landscape mean.
  - Knots are fixed at fit time.
  - AIC and AICc at nominal ν over-select flexibility.
  - Describe the monotonicity report.
  - Tails outside the sampled range are poorly identified.
- **`?terradish_scale_optim`:** units; the `"terradish"` default versus multiScaleR; σ counted in K; the search is secondary to the joint model.

### Measurement models

- **`?mlpe`**
  - It is ML, not REML, so AIC comparisons across conductance formulas are valid.
  - ρ is on the logit scale, with the boundary handling from I6.
  - In package validation, model-based SEs understated the spread among replicates about 1.7-fold.
- **`?mlpe_covariates`, `?pairwise_endpoint_covariates`, `?pairwise_covariates`**
  - Surfaces are conditional on the pairwise terms, and signs can change with the mean structure; report fits with and without them.
  - A `~ 1` model is IBD plus IBE, not IBE alone.
  - Document the new input function, and the IBE:IBR ratio (I12).
- **`?wishart_covariates`**
  - The kernel K_k = −½ H D_k H, built from the same pairwise transforms as `mlpe_covariates()`.
  - The distance-form identity with the MLPE mean (F8b), and the IBE:IBR ratio (I12).
  - The contrast form for covariance input; λ_k ≥ 0 versus the unconstrained γ_k; chi-bar-square tests and boundary reporting (I11).
  - **Interpretation.** A positive λ_k means genetic differentiation increases with environmental difference, conditional on resistance. It is not specific to one mechanism: assortative mating, dispersal filtering and local adaptation all produce it, and so can site-level variance that rises at environmental extremes. The same is true of the MLPE pairwise term.
  - At ν equal to the SNP count, AIC selected the kernel in 19 of 40 null simulations, and in none at the replicate-calibrated ν. Recommend `terradish_rescale_nu()` sensitivity and fixed-nuisance CV (C10).
- **`?generalized_wishart`, `?check_distance_response`**
  - Use admissible, covariance-derived squared distances only.
  - FST ratio estimators have no Wishart ν.
  - Include the ν paragraph.
- **`?wishart_covariance`**
  - The contrast likelihood, and its equivalence to GW on `dist_from_cov(C)`.
  - Invariance to centering.
  - The full-rank form was removed, and why.
  - The diagonal-mode caveat and the ν paragraph.

### Inference and model comparison

- **`?terradish_rescale_nu`:** what changes (SEs, tests, AIC) and what does not (estimates). Report conclusions across a plausible range of ν.
- **`?terradish_ibe_ratio`:** the definition, its equivalence across families, and the delta-method interval.
- **`?aic_table` and the anova methods**
  - The comparison contract, the new columns, the nesting rules, and the `nobs` and AICc conventions (I1).
  - Information criteria at nominal ν over-select extensions. Use CV to select conductance terms, and `terradish_rescale_nu` for sensitivity.
- **`?terradish_cv_folds`, `?terradish_folds`:** the interpretation block from Section 8, per-formula models, the two `nuisance` modes (C10), paired differences, the checkpoint signature, and failed folds.
- **anova methods:** the chi-bar-square reference for λ (I11).

### Design, utilities and re-exports

- **`?simulate_covariance_response`, `?covariance_response_power`:** map-unit σ and `nu_fit`. Results are internal to the model, inherit any error in ν, and cannot reveal that error.
- **`?terradish_assess_settings`:** a speed probe, not a convergence or statistical check.
- **`?slim_terradish`:** `vcov` and `confint` now work; slim fits cannot predict.
- **`?scale_covariates`, `?scale_to_0_1`:** the `reference` argument, and per-layer behavior.
- **`?crop_to_focal_buffer`:** the caveat and a sensitivity recipe.
- **Optimizer controls:** the new stopping rule.
- **Re-exports (`R/reexports.R` header):** point to landgraph, and carry the caveats into the vignettes.

### Acceptance for Phase 6
- `R CMD check` shows no Rd warnings.
- Examples pass under `--run-donttest`.
- `inst/WORDLIST` is updated and `spelling::spell_check_package()` is clean.

---

## 11. Phase 7: vignettes, README, NEWS, DESCRIPTION (commit C7)

Standards for every vignette:
- explain how to interpret the outputs;
- annotate the code;
- end with a quick reference and links to the other vignettes;
- re-knit from scratch after Phases 2–5 (delete `vignettes/*_cache`).

### getting-started.Rmd
- "What is terradish for?" and "A conceptual road map": the 0.1.0 scope.
- The new default, `directions = 8`.
- "Choosing a measurement model": MLPE and least squares versus Wishart; covariance and distance inputs give one contrast likelihood.
- A new subsection, "What ν means" (Appendix C).
- Coefficients: relative conductance, plus `vcov()` and `confint()` (line 864 already advertises them).
- "Support-constrained prediction" (lines 737–770): re-run after F3, and state exactly what focal support is.
- A new subsection, "Convergence and starting values" (I2–I4).
- Update the quick reference.

### gaussian-scale-optimization.Rmd
- "What sigma means here" and its limits (lines 47–96): what the bounds mean, dispersal averaging, β per SD of the smoothed layer, the truncated kernel and the default upper bound.
- Steps 3–7: re-run. Every number changes, because the melip rasters have even dimensions.
- Step 4: the new summary and the profile interval.
- Step 5: add a `terradish_cv_folds` comparison (C1).
- Update "Important current limitations".
- Add a short outer-search section on units and kernels.

### spline-conductance.Rmd
- Basis section (lines 54–121): unpenalized, fixed df, no partition-of-unity claim, centered basis.
- "Model comparison" (line 348) and the AICc stepping advice (around line 786): replace with CV (C1) and the monotonicity report.
  - Explain that AIC over-selects.
  - Explain that a selected spline can reflect how the process maps onto the graph rather than an ecological effect.
- Tail clamping (lines 447–472): re-run after F3 and F4.

### model-comparison.Rmd
Preserve the owner's P0.1 text.
- LRT section: nesting rules (I8) and the ν caveat.
- Information criteria (lines 225–314): over-selection at nominal ν, and `terradish_rescale_nu` sensitivity.
- Write the new CV section where S10 left the comment: `terradish_folds` and `terradish_cv_folds`, per-formula models, paired differences, and repeats (if C7 is done).
- "Combining the evidence":
  - use CV to select conductance terms;
  - use ν sensitivity for inference;
  - treat LRT at a fixed ν as descriptive only.
- Update the quick reference.

### wishart-covariance.Rmd
- Overview: the contrast likelihood, and why the full-rank form was removed.
- Preparing data:
  - diagonal modes (`gower` versus `within`);
  - rare-variant weighting, and a recommendation to filter by MAF;
  - unequal n means unequal sampling variance.
- "Choosing nu" (lines 159–238): rewrite.
  - ν is set by the user.
  - Linkage and shared history reduce the information.
  - The block jackknife gives an upper bound.
  - Use CV and rescaling.
- "The role of nu" (lines 364–441): tie it to `terradish_rescale_nu`.
- Fix the caption at line 293.
- GW section: the FST caveat.
- A short section on environmental-difference terms that points to the IBE vignette.

### ibe-ibr-workflow.Rmd
Rebuild it around one `pairwise_covariates()` object used by both families:
- geographic distance and scaled environmental differences in one object;
- fit IBR + IBE with `mlpe_covariates()` (MLPE) and with `wishart_covariates()` (Wishart, on a covariance-derived response). Report the IBE:IBR ratio from both side by side, and explain why the ratios are comparable while the raw coefficients are not;
- relabel "IBE only" as "IBD + IBE (uniform surface)";
- the conditional-surface caveat, with θ shown both with and without the pairwise terms;
- model comparison: LRT (χ² for γ, chi-bar-square for λ), `terradish_rescale_nu()` sensitivity for the Wishart fit, and fixed-nuisance CV of IBR against IBR + IBE (C10);
- interpretation limits: the term is not mechanism-specific, and λ ≥ 0 while γ is unconstrained;
- keep the visualization.

### large-landscapes.Rmd
- Solvers: `direct`, `auto` and `amg`; the `auto` rule; AMG's sensitivity to contrast.
- **Timing table** (around lines 115–123): delete it. Keep it only if the owner generates it with a committed `inst/benchmarks/large-landscape-scaling.R`, stores the raw outputs in `inst/benchmarks/results/`, and reconciles it with DESCRIPTION.
- A new section, "Cropping and coarse warm starts": the R5 numbers, and exact refinement as a requirement.
- Mark worker reuse as experimental (R4).

### simulation-design.Rmd
- In the ν power section (lines 197–288), add `nu_fit` and the caveat that results are internal to the model.
- Regenerate the power `.rds` files with `precompute.R` after F6 and R9.
- **Build-ignore problem:** `.Rbuildignore` excludes `vignettes/*.rds`, so CRAN builds render these sections empty. Either ship compact CSV summaries in `inst/extdata/` and read those, or drop the precomputed sections.

### README
- Rewrite it around Section 1.1.
- Fix the `pairwise_endpoint_covariates` example (line 334).
- Add sections titled "Interpreting estimates", "Experimental features" (with branch install instructions) and "Changes in 0.1.0".
- Keep the DOI.

### NEWS.md (0.1.0)
- **Scope:** the moved features, with reasons.
- **Breaking changes:**
  - `wishart_covariance` now uses the contrast likelihood.
  - `wishart_covariates` is reformulated: contrast form for covariance input, kernels from the MLPE pairwise transforms (default `absdiff`), `normalize` removed, pairwise objects accepted.
  - The default is `directions = 8`.
  - The solver set is reduced.
  - `exact_refine = FALSE` is removed, and `terradish_grid` accepts only `"none"`.
  - The legacy CV functions are removed.
  - `terradish_scale_optim` has a new default kernel and counts σ in df.
  - The simulator takes σ in map units.
  - The Gaussian `sigma_upper` default changes.
  - Spline `conductance()` values shift by a constant (F4).
- **Bug fixes:** F1, F3, F7, I3, I5, I6, I8, I9, I11, C2–C4, R2, R3, R8.
- **New:** N1–N15 as implemented, and `scale_covariates(reference = )`.
- **Documentation:** the ν guidance, the vignette rewrites, and the removal of out-of-core claims and unreproduced benchmarks.

### DESCRIPTION and metadata
- **Version and date:** 0.1.0 and the release date.
- **Title:** keep it, or change it to "Likelihood-Based Estimation of Landscape Conductance Surfaces from Genetic Data" (**owner decision**).
- **Description:** use the text in Appendix C.
  - Remove the diagonal-variance extension and the selected-pair analysis. Describe the environmental terms as shared by both families.
  - Remove "fewer than a million cells" unless a committed benchmark supports it.
- **Imports:** re-verify each one after the removals (`ggplot2`, `MASS`, `multiScaleR`, `nlme`, `parallel`, `splines`, `utils`, `grDevices`, `methods`).
- **landgraph:** require the minimum version from Section 14.
- **Other files:**
  - Update `CITATION.cff`, `.zenodo.json`, `cran-comments.md` (first release; the status of the landgraph dependency) and `inst/WORDLIST`.
  - Check `inst/benchmarks/compare-radish-terradish.R`.
  - Re-run `inst/examples/*.R`.

### Acceptance for Phase 7
- All vignettes build.
- The Phase 1b grep is still clean.
- No stale numbers remain.

---

## 12. Phase 8: validation and release

1. **Checks.**
   - Run `devtools::document()`, `devtools::test()`, `devtools::check(args = "--as-cran")`, and a final vignette-enabled check in RStudio.
   - Target: 0 errors, 0 warnings, and only the NOTEs for a new submission and for landgraph.
2. **Numerical report in `log.md`.** Rerun `dev/release_0.1.0/core_baseline.R` with the same explicit settings.
   - Unchanged: (a), (b), (c), (f) and (h).
   - Changed by design: (d) and (g).
   - Same θ and logLik: (e).
   - Any other difference blocks the release until it is explained.
3. **Audit scripts as regression tests.** Fix their paths, then rerun:
   - `core/exp1b_wishart_clustered.R`, `core/exp7_optim.R`, `core/exp9_rho0.R`;
   - `scale/shift_check.R`, `scale/ksr_units.R`;
   - the spline-support experiment in `flex/`;
   - `cv/e3.R` (using the `terradish_cv_folds` equivalent), `cv/e5.R`, `cv/e6.R`;
   - the `compute/` checks for `scale_to_0_1` and the crop example;
   - `ext/check1`, `ext/check2` and `ext/check5` (rerun check5 with the new `absdiff` default and report the result), and `cv/e2.R` under both `nuisance` modes.

   Each must show the corrected behavior. Save the outputs.
4. **Coverage.** Run `covr::package_coverage()`. Aim for at least 80% coverage on the files touched in Phases 2–5.
5. **CRAN tests.** After the split, four test files skip entirely on CRAN. Keep a fast subset (under 2 s) of the derivative and AMG correctness tests running on CRAN.
6. **Release.** After owner sign-off:
   - set Version 0.1.0;
   - merge `release/0.1.0` into `master`;
   - tag `v0.1.0`;
   - run the Section 2 sync.

---

## 13. `experimental` branch (Phase 9)

`EXPERIMENTAL.md` lists the moved features, their status, and these known issues.

**Directed model**
- With `edge_gradient` covariates the chain is reversible for every γ, so only θ/2 − sγ is identified.
- `edge_flow` is untested.
- Needed: an identifiability analysis for non-gradient fields, and a guard against the same covariate entering both terms.

**Hierarchical field**
- `conductance_model` is ignored (`R/hierarchical_conductance.R:259,295,536`).
- The full Hessian is formed at every τ² grid point.
- It uses a conditional AIC with a naive edf.
- `aic_table` rejects these fits.
- Non-convergence is silent.
- There is no restricted spatial regression.
- Under the correct model and ν it is not over-selected, but it absorbs process misspecification.
- Adding it inflated the covariate SE about 5-fold (`package_audit/identifiability/field_info.R`).

**Drift covariates**
- A site-specific diagonal is algebraically identical to site-specific additive distance. A local halo of low conductance, sampling variance and within-deme diversity therefore all load on it. In `identifiability/halo.R`, a site constant explains 99.9% of a halo's effect.
- The sign reverses under rare-variant weighting.
- There is a 1/n_i artifact.
- Needed: sampling-variance offsets, a MAF filter, and a new name.

**Pair subsets**
- Wrapping `mlpe_covariates` drops γ.
- There is no CV subsetter.
- The order guard is too strict.

**Other items**
- **`linear_conductance`:** constrain one coefficient, or the scale.
- **Kron reduction:** not out-of-core, and its covariance output is not the fitting kernel.
- **`block_cg`:** fails with six or more right-hand sides.
- **Legacy CV:** external-θ scoring (`R/radish_cv.R:691`), and φ re-profiling that favors measurement extensions.
- **`terradish_grid` screening approximations:** there is no exact refinement.

---

## 14. landgraph companion changes

Repository: `C:\Users\peterman.73\OneDrive - The Ohio State University\R\Packages\landgraph`.

**L1. Attributes on covariance output.**
- `cov_from_biallelic` and `cov_from_genetic_data` record `attr(x, "diagonal")`, and set `attr(x, "centered") = "sites"`.
- This is not a breaking change.
- *Test:* the attributes are present.

**L2. Documentation.**
- The within-population diagonal is on a different scale and is unsuitable for Wishart likelihoods.
- Per-locus standardization up-weights rare variants on the diagonal; recommend a MAF filter.
- Unequal n means unequal sampling variance.
- `fst_from_biallelic` output is not a Wishart response.
- `pca_dist` fills missing data with the mean.

**L3. Decision gate.**
- Run the audit's `core/exp11_within.R` design with the new contrast `wishart_covariance`.
- If `diagonal = "within"` still biases θ (a wrong sign, or bias above 2 SE), the **owner decides**: either a `"gower"` default for grouped data, or keep the current default plus the terradish warning (R12).
- Record the result and save the output.

**L4. Versioning.**
- Bump the landgraph version, and have terradish require it.
- A CRAN release of terradish requires landgraph to be on CRAN first.

`edge_gradient` and `edge_flow` stay in landgraph. terradish 0.1.0 no longer re-exports them.

---

## 15. Downstream work for the owner (not for the package worker)

After `v0.1.0` is released:
- Refit SLiM experiments E1, E2/E2N, E3 and E8 with the release.
  - The E2 truths used the old crop, but a one-cell shift is negligible at σ = 18.
  - Re-estimate the E2N σ̂.
- Remove E4–E6 and their nulls from the paper.
- Keep E7/E7N and refit them with the reformulated kernel. The SLiM mechanism weights mating and dispersal by Δe², so use `transform = "sqdiff"`, or report both transforms. Add the IBE:IBR ratio.
- Refit the Sceloporus Gaussian models.
- Replace AIC at nominal ν with CV, plus sensitivity to ν via `terradish_rescale_nu`.
- Update the function table, Table 5, Figures 3–4 and the supplement.

---

## Appendix A. File checklists

**Deleted on `master` (C1b):**
- `R/`: `directed_conductance.R`, `hierarchical_conductance.R`, `wishart_drift_covariates.R`, `pair_subset_measurement.R`, `radish_multiscale.R`, `kron_reduction.R`
- `src/`: `directed_sparse_lu.cpp`
- `tests/testthat/`: `test-directed-conductance.R`, `test-hierarchical-conductance.R`, `test-wishart-drift-covariates.R`, `test-pair-subset-measurement.R`, `test-kron-reduction.R`, `test-block-cg.R`, `test-cv-helpers.R`
- `vignettes/`: `directional-conductance.Rmd`, `hierarchical-conductance.Rmd`, `vignette-directional.rds`, `vignette-hierarchical.rds`
- `dev/`: `make_vignette_directional.R`, `make_vignette_hier.R`

**Created in C1a:** `R/aic_table.R` and `R/cv_helpers.R` (and `R/utils_matrix.R` if A4 uses it). The helpers relocated into `R/radish_grid.R` and `R/terradish-package.R` are listed in Section 4.

**Edited for removal hooks (C1b):**
- `R/`:
  - `radish_cv.R` (or deleted) and `aic_table.R`;
  - `radish_optimize.R`, `radish_algorithm.R`, `radish_conductance_model.R`, `radish_grid.R`;
  - `mlpe_covariates.R` (links only), `generalized_wishart.R`, `wishart_covariance.R`, `covariance_response_power.R`;
  - `terradish-package.R`, `reexports.R`, `RcppExports.R`.
- `src/RcppExports.cpp`
- Tests: `test-analytic-derivatives.R`, `test-verbose-and-tau2.R`, `test-model-selection.R`, `test-ibe-workflow.R`, `test-landmark-grid-cv.R`, `test-spatial-cv.R`, `test-legacy-wrappers.R`
- Vignettes: all eight kept vignettes (links and removed sections) and `precompute.R`
- `README.md` and `inst/benchmarks/synthetic-solver-scaling.R`

**New test files by the release:**
- `test-gaussian-alignment.R`, `test-spline-support.R`
- `test-inference-methods.R`, `test-optimizer-stopping.R`, `test-nu-rescale.R`, `test-anova-nesting.R`
- `test-cv-folds-core.R`
- `test-prediction-new-surface.R`, `test-pairwise-covariates.R`, `test-power-nu-fit.R`
- Rewritten `test-wishart-covariance.R` and `test-wishart-covariates.R`
- `test-ibe-ratio.R`, `test-cv-fixed-nuisance.R`
- Updates to `test-scale-optim.R`, `test-gaussian-scale-conductance.R`, `test-smooth-loglinear-conductance.R`, `test-terra-workflow.R` and `test-covariance-response-power.R`

## Appendix B. Audit evidence index (`review_work_20260924/package_audit/`)

| Topic | Scripts | Plan items |
|---|---|---|
| wishart_covariance collapse; diagonal modes | `core/exp1_wishart.R`, `core/exp1b_wishart_clustered.R`, `core/exp1c_biallelic.R`, `core/exp10_covgen.R`, `core/exp11_within.R` | F2, R12, L3 |
| MLPE is ML; ρ at the boundary | `core/exp2_mlpe.R`, `core/exp9_rho0.R` | I6 |
| ν scaling; missing vcov | `core/exp3_s3.R` | I1, I7 |
| linear_conductance | `core/exp4_linear.R` | S5 |
| Focal-support spline bug | `core/exp5_support.R`, `flex/` | F3 |
| Boundary starts | `core/exp7_optim.R` | I4 |
| Gaussian crop; σ bounds and units | `scale/shift_check.R`, `scale/smoother_check.R`, `scale/bound_and_units.R`, `scale/ksr_units.R`, `scale/exp_outer.R` | F1, F5, F6, I9 |
| Kernel, drift, pairwise covariates, anova | `ext/check1`–`check6` | S3, S6, F8, I8, I11, I12, R8 |
| CV; ν invariance; power coverage | `cv/e1.R`–`cv/e7.R` | C1–C8, C10, R9 |
| Solvers, crop bias, scale_to_0_1, parallel, slim | `compute/e1`–`e8` | R2–R5, F7, S8, S9 |
| Identifiability | `identifiability/field_info.R`, `identifiability/halo.R` | Section 13 |

## Appendix C. Canonical text

**Standard ν paragraph:**

> For Wishart likelihoods, `nu` scales the log likelihood: point estimates do not depend on it, but standard errors, likelihood-ratio tests, and information criteria do. `nu` is a dispersion parameter that the user sets, not the number of markers. Linkage and, more importantly, the population history that all loci share make the effective information far smaller than the SNP count; in forward-time simulations used to validate terradish, standard errors at `nu` equal to the SNP count were three to five times smaller than the spread among replicate histories. Treat model-based standard errors as lower bounds, use `terradish_cv_folds()` to choose among conductance formulas (its rankings do not depend on `nu`), and use `terradish_rescale_nu()` to report how inferences change across plausible values.

**DESCRIPTION `Description:` (draft):**

> Maximum likelihood estimation of landscape conductance surfaces from genetic data on raster graphs, using Newton and quasi-Newton optimization with analytic derivatives. Fits log-linear conductance models with optional joint estimation of Gaussian smoothing scales and spline response shapes, linked to pairwise genetic distances through maximum likelihood population effects models or to allele-frequency covariance through a Wishart likelihood on site contrasts, with the same pairwise distance and environment covariates available in both. Provides exact sparse and algebraic multigrid solvers for large rasters, fixed-domain spatial cross-validation for comparing conductance formulas, and simulation tools for study design. Works natively with 'terra' objects.
