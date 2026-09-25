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

## C1b: core split, version 0.0.50

Removed experimental features and their exclusive documentation/tests while
preserving shared numerical helpers and internal compiled code permitted by S9.
All nine frozen baseline results are identical. Scope and globals checks pass.
Source tests: 703 passed, zero failures/errors, one expected legacy warning
and one intentional parallel skip. Full vignette-enabled package check:
Status OK, zero errors, warnings, or notes. See phase1b_00check.log.
Metadata-only version bump followed validation.

## Phase 3 acceptance investigation

The initial outer-optimizer prototype reports code 2 for the corrected even-grid
Gaussian baseline, projected gradient 0.0007226463. The other eight fits report
code 0. The owner authorized investigating and fixing the cause while retaining
the specified stopping thresholds. This prototype does not yet exercise the
new inner warm-start path; full-package validation remains required.

## C2: correctness fixes, version 0.0.51

Implemented F1-F8 and dependency L4. The complete candidate passed 803 assertions
with zero failures/errors, 10 intentional Gaussian-support warnings, one legacy
maxit warning, and one intentional parallel skip. Full vignette-enabled check
completed Status OK: zero errors, warnings, or notes. See phase2_00check.log.
Actual-candidate baseline comparison passed. MLPE, least squares, generalized
Wishart, odd-grid Gaussian, and MLPE covariate results are unchanged. Spline
coefficients/SE/likelihood are unchanged; nuisance scaling absorbs centering.
Covariance Wishart now matches generalized Wishart. Even-grid Gaussian changes
are intentional. Kernel lambda rescales by the removed kernel variance.

Landgraph L3 audit and the owner-approved gower default are complete. Companion
commit 36a4e17 (0.0.3) passed 74 assertions and full check. Installed only in
an isolated validation library. Ordinary library remains unchanged. No push
or merge occurred. The dependency must be released before CRAN submission.

The Phase 2 runner saved its test/check RDS receipts before a trailing top-level
else parse error in reporting. Both complete R CMD check and structured test
results independently establish the results above. Runner corrected after both
processes finished; no running script was edited.

## C3-C7: completed implementation and documentation

C3 (.52) adds inference methods and corrects numerical stopping, with the
owner-authorized Gaussian repair preserving the specified thresholds. C4 (.53)
adds fixed-domain predictive CV, common-fold accounting, and checkpoint safety.
Their detailed evidence is in `phase3_phase4_validation.md`. The exact C4
archive check passed with zero errors, warnings, or notes.

C5 (.54, 9b75ff4) implements robustness and prediction. C6 (.55, e81acdb)
revises function documentation. Both checks are clean; C6 skips redundant tests
on unchanged code. See `phase5_phase6_validation.md`. C7 (.56, ab6c020) rewrites
README and all eight core guides. All guides and standalone examples run; the
examples/vignette check has zero errors, warnings, or notes. The small power
example yields no defensible sample-size recommendation. See
`phase7_validation.md` for the precise scope of verification.

## Final numerical regression report (C8 candidate)

The final nine fits use the original explicit four-direction baseline settings.
All have convergence code zero. `core_phase8.rds`, `phase8_baseline.log`, and
`phase8_baseline_comparison.csv` retain parameters, nuisance estimates, standard
errors, and comparisons with both .47 and the corrected C2 reference.

| Model | Original log likelihood | Final log likelihood | Interpretation |
|---|---:|---:|---|
| a: MLPE | 2160.184078 | 2160.184078 | Unchanged |
| b: least squares | 1960.571248 | 1960.571248 | Unchanged |
| c: generalized Wishart | -609.884981 | -609.884981 | Unchanged |
| d: covariance Wishart | -664.826513 | -609.884981 | Intended site-contrast correction |
| e: spline | 2182.350172 | 2182.350172 | Same theta and likelihood; centered basis |
| f: Gaussian, odd | 1517.628227 | 1517.628227 | Unchanged within tolerance |
| g: Gaussian, even | 1549.991372 | 1517.537750 | Intended alignment correction |
| h: MLPE covariates | 2222.044115 | 2222.044115 | Unchanged |
| i: Wishart covariates | -609.880167 | -609.880167 | Historical sqdiff comparison; lambda rescales |

The largest final-minus-C2 differences are 9.73e-10 in log likelihood,
2.85e-6 in theta, 8.33e-6 in nuisance parameters, and 5.66e-6 in standard
errors. All satisfy the declared tolerances. The kernel baseline deliberately
retains sqdiff to compare with the old kernel; separate default-absdiff audits
are in `phase8_audit_interpretation.md`. No unexplained numerical difference
was found.

The final noise-CV audit gives 3/8 substantive Wishart and 2/8 MLPE selections
under fixed nuisance scoring, versus 8/8 for both when nuisance parameters are
reprofiled. One Wishart difference is a numerical tie. Failed model and
baseline folds remain explicit. The detailed audit report also records the
remaining AIC and near-zero-rho limitations.

The experimental branch now has documentation commit 38f57d8 (.48), with its
implementation unchanged from .47. Neither branch has been pushed or merged.
The landgraph CRAN prerequisite and owner release approval remain open.

## C8 final local validation and adoption (0.0.57)

Adopted phase8-final-package after its caption-corrected archive check passed
0 errors, 0 warnings, 0 notes (12m21.7s). Tests were deliberately skipped in
that final documentation check; byte-identical R/C++/test sources passed the
complete .57 check (27m29.3s), including 990 installed assertions, 197 diagnostic
warnings, and no skips. Clean source tests passed 988 assertions with the same
warnings and one intentional parallel skip. All nine numerical baselines meet
C2 tolerances. Eight fast CRAN assertions passed in 1.21 seconds combined.

R line coverage: 86.79266% overall, 88.02621% across 25 changed files. Individual
coverage below 80% is disclosed in RELEASE_REVIEW.md. Native gcov collection
failed with exit code 6; no native line-coverage result is claimed. Earlier
interrupted attempts are distinguished in phase8_validation_environment.md.
The final archive contains all 43 figures across eight guides. Representative
plots were inspected; no complete browser-layout inspection is claimed.

Final archive SHA-256:
dade1f8982f39ae55cf9e5f8a5b32cc87e0c3d55c5140cb74801c4797981d5cc.
See phase8_final_fingerprint.json. No numerical implementation changed in C8.
The interactive RStudio check (or owner acceptance of its command-line
equivalent), companion landgraph publication before CRAN submission, and owner
release sign-off remain open. No master merge, push, final release tag, or
submission has been performed.

Final diff cleanup corrected two obsolete setup comments in the spline and
model-comparison vignette sources. This comment-only change followed the final
archive check; evaluation settings, executable chunks, and rendered prose are
unchanged. The archive fingerprint refers to the preserved checked candidate.

## September 25: documentation math correction (0.0.58)

Repaired encoding damage, equation notation, and symbol/column consistency.
Preserved the owner's existing model-comparison paragraph wrapping. All eight
vignettes rebuilt; 62 Rd checks passed; 61 vignette and 124 help math expressions
rendered without errors. Executable R syntax is unchanged. No complete numerical
suite or archive check was repeated. See ../docs_math_20260925/REVIEW.md.
