# terradish 0.1.0 implementation log

## Phase 0: 2026-09-24

Source is `6a2c13f` (0.0.47). The original vignette edit only wraps existing
paragraphs. The owner approved including it in C0. Created `release/0.1.0`,
`experimental`, and tag `v0.0.47`. The experimental branch retains all original
features; its explanatory document is still pending.

R 4.6.1 and installed landgraph 0.0.1 were verified live. Baseline settings:
four-neighbor graph, exact curvature, direct solver, fixed simulation seed
20260924, and 200 permitted Newton iterations. All 37 melip sites are used on
the original rasters for models a-e and h-i. Gaussian models use a six-cell
aggregate and complete-case sites; model f crops the aggregate to odd raster
dimensions. Explicit scale bounds satisfy the plan's kernel-support rule.

All nine models completed without warnings or errors. Full values, iteration
counts, nuisance estimates, coefficients and SEs are in the frozen RDS file;
the console record is `baseline.log`.

| Model | log likelihood | Maximum absolute gradient | Seconds |
|---|---:|---:|---:|
| a: MLPE | 2160.1841 | 5.15e-9 | 5.35 |
| b: least squares | 1960.5712 | 7.48e-11 | 2.78 |
| c: generalized Wishart | -609.8850 | 1.40e-10 | 2.71 |
| d: covariance Wishart | -664.8265 | 3.01e-9 | 2.26 |
| e: spline | 2182.3502 | 9.74e-12 | 7.92 |
| f: Gaussian, odd | 1517.6282 | 13.08031 | 4.07 |
| g: Gaussian, even | 1549.9914 | 1.72e-8 | 3.31 |
| h: MLPE covariates | 2222.0441 | 2.11e-11 | 9.46 |
| i: Wishart covariates | -609.8802 | 8.34e-11 | 3.62 |

Model f's scale is at its lower bound (0.135 map units). Its large raw gradient
is in that constrained coordinate; the largest remaining component is
1.014345e-5. Preserve this evidence when validating projected-gradient stopping
in Phase 3. None of these fits is on the no-structure measurement boundary.

The base namespace globals scan returned no unresolved function names.

Base source tests: 920 passed, zero failures or errors, one expected legacy
wrapper warning and one intentional parallel-test skip. Elapsed time was
602 seconds. See `base_test.rds` for structured results. The runner later
reported a parse error because its source file was revised while R was still
reading it; the test results had already been saved. The saved runner now
parses cleanly. Do not edit a running validation script.
Base vignette-enabled package check completed with Status: OK (zero errors, warnings or notes); see `base_check.log` and `base_00check.log`. The check took 52m42s. After check completion, the early runner revision caused a trailing parse error before its RDS receipt was saved; the complete R CMD check log independently records Status: OK.
The first build attempt failed while copying deep `.git/refs/codex` paths on
Windows, before package checks. Its log and RDS are retained with the suffix
`copy_failure`. The replacement builds a temporary copy of working-tree
package files, without `.git`, and still builds the vignettes.

## Companion landgraph work

Owner approved correcting L1's inaccurate unconditional `centered = "sites"`
label. Live examples produced nonzero row sums for both unequal-sample
biallelic covariance and within-diagonal covariance. Added construction labels
without changing values: `sites`, `pooled_allele_frequency`, and
`sites_before_diagonal_replacement`. Biallelic diagonal construction is labeled
`normalized_dosage`; genetic-data covariance retains `gower` or `within`.

L1 and L2 covariance/FST documentation changes are committed as `501aa23` on
`landgraph` branch `codex/terradish-core-companion`, version 0.0.2. All 73 tests
passed; vignette-enabled `--as-cran --no-manual` check returned Status: OK.
The check preceded the metadata-only version bump, as required by the phase
workflow. `pca_dist()` is in terradish, where mean imputation is already
documented. L3's contrast-model experiment and default decision remain pending.
No package was installed, pushed, or merged; base terradish validation still
uses installed landgraph 0.0.1.

## C1a: helper relocation, version 0.0.49

All 414 function bodies, formals, and classes are unchanged. All nine frozen
numerical baselines are identical. Source tests: 920 passed, zero failures or
errors, one expected legacy warning and one intentional parallel skip.
The full vignette-enabled check completed with Status: OK, zero errors,
warnings, or notes, in 53m48s. See phase1a_00check.log and phase1a_check.rds.
The check used the frozen source before the metadata-only version bump.
