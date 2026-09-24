# terradish core overhaul: owner review

The core implementation and documentation follow the supplied 0.1.0 plan and
the four decisions recorded during this task. Implementation and local
validation are complete. The release gates below remain open.

## What changed

The core now fits symmetric resistance models with log-linear, spline, and
joint Gaussian conductance factories. Covariance Wishart uses site contrasts.
MLPE and Wishart environmental extensions accept the same pairwise inputs,
with absolute differences as the default endpoint transform. Gaussian grid
alignment, spline support, scaling, nested comparisons, and optimizer stopping
were corrected. New inference methods, convergence records, predictive CV,
nu sensitivity, retained transformations for new rasters, and power-study
nu_fit controls are documented and tested.

README and all eight retained guides were rebuilt around this scope. Removed
research features remain on `experimental`, with their known limitations in
EXPERIMENTAL.md. That branch preserves the pre-overhaul implementation; it
does not receive the core repairs automatically.

The companion landgraph now records truthful covariance-construction metadata,
documents sampling and diagonal limitations, and uses the owner-approved Gower
default for grouped data. Its numerical output was not changed merely to add
metadata.

## Local commit record

| Batch | Version | Commit |
|---|---|---|
| Preparation, preserving the owner's vignette wrapping | .48 | c80425d |
| Helper relocation | .49 | 488ac6d |
| Core split | .50 | 604a3a1 |
| Correctness repairs | .51 | c3ac708 |
| Inference and numerical stopping | .52 | c2aa209 |
| Cross-validation | .53 | 0bfb63e |
| Robustness and prediction | .54 | 9b75ff4 |
| Function documentation | .55 | e81acdb |
| Guides, examples, and release metadata | .56 | ab6c020 |
| Final validation, fast CRAN tests, and complete guide outputs | .57 | This commit |
| Experimental-branch notes | .48 | 38f57d8 |
| landgraph metadata/documentation | .2 | 501aa23 |
| landgraph grouped-data default | .3 | 36a4e17 |

No push, master merge, release tag, or CRAN submission has been performed.
The original master commit and v0.0.47 checkpoint remain intact.

## Validation evidence

- C4's exact archive and C5's complete package check pass with zero errors,
  warnings, or notes. C6 and C7 examples/vignette checks also pass cleanly;
  their redundant test runs were skipped because executable R code was unchanged.
- All eight revised guides render. README and IBE example fits converge. The
  standalone examples complete, including the power example; its fit rates
  are 0.8-1 and no small design reaches its recovery target.
- All nine final numerical baselines converge and meet the declared tolerance
  against the corrected C2 reference. The covariance and even-grid Gaussian
  changes are intentional. See log.md and phase8_baseline_comparison.csv.
- Independent contrast-Wishart and dense MLPE calculations agree with the
  implementation within 1e-8. Alignment, support, scaling, crop, checkpoint,
  internal-scale CV, nesting, and boundary diagnostics are retained in the
  final audit reports.
- Eight fast derivative/AMG assertions run on CRAN and passed in 1.21 seconds
  combined, excluding package loading/compilation.
- Final source suite: 988 assertions passed, zero failures/errors, 197 diagnostic
  warnings, one intentional parallel skip. Installed archive tests passed 990
  assertions with the same warnings and zero skips, including worker consistency.
  The .57 archive check passed with zero errors, warnings, or notes, including
  installed tests, examples, donttest examples, and all eight vignette rebuilds
  (27m29s for the check stage).

Inspection of the packaged plots found that two old global NOT_CRAN switches
suppressed spline/model-comparison guide outputs in ordinary builds. Explicit
renders and check rebuilds had evaluated them. Those switches are now removed;
the final archive contains all 43 guide figures. A stale spline caption was
also corrected to describe the actual axes and descriptive regression line.
The final caption-corrected archive check passed with zero errors, warnings,
or notes in 12m21.7s, including examples and all eight vignette rebuilds.
Tests were skipped in that check because its R, C++, and test sources are
byte-identical to those in the complete passing .57 check. Representative
covariance, IBE, Gaussian-kernel, spline, likelihood-surface, and conductance-map
images were inspected directly. This is plot inspection, not a claim that all
43 figures or the full browser layout were visually inspected.

R-source line coverage is 86.79% overall and 88.03% across the 25 R files
changed in Phases 2-5. The aggregate exceeds the 80% aim, but three individual
changed files fall below it:

| File | Covered lines | Coverage |
|---|---:|---:|
| R/generalized_wishart.R | 106/166 | 63.86% |
| R/mlpe.R | 160/209 | 76.56% |
| R/radish_graph.R | 155/198 | 78.28% |

The full test suite ran under R tracing with normal native compilation.
Detailed gaps are saved in phase8_r_uncovered_touched_lines.csv. The separate
native-coverage attempt failed during gcov collection with exit code 6;
native line coverage is unverified. This collector error does not invalidate
the independently passing numerical tests.

Source-test warnings comprise 151 unused-covariate notices, 22 iteration-limit
notices, 11 Gaussian support warnings, six model-ranking/boundary notices, four
optimizer stalls, two duplicate-cell warnings, and one failed uniform-baseline
score. Their presence is not summarized as a warning-free source test run.

The validation environment report distinguishes completed checks from attempts
interrupted for thread contention and shared instrumented build artifacts.
No interrupted attempt is counted as a pass.
The final devtools check used its no-manual default and disabled remote
incoming checks. This establishes local validation, not CRAN acceptance.
Final archive: dev/check/phase8-final/terradish_0.0.57.tar.gz. Its SHA-256 is
`dade1f8982f39ae55cf9e5f8a5b32cc87e0c3d55c5140cb74801c4797981d5cc`;
candidate source fingerprints are in phase8_final_fingerprint.json.
The initial complete check and its archive remain preserved separately under
phase8_check.* and dev/check/phase8.

## Interpretation limits retained in the release

Fixed predictive scoring selects a noise environmental extension in 3/8
Wishart and 2/8 MLPE replicates; reprofiled scoring selects it in 8/8 for both.
One Wishart comparison is a numerical tie. These eight-replicate results are
regression evidence, not calibrated false-positive rates. AIC can still favor
an environmental term under misspecification. Environmental coefficients and
IBE:IBR ratios remain conditional associations, not identified mechanisms.

Near rho = 0, the MLPE audit matches least squares and the outer fit converges,
but an inner nuisance stall is reported. Boundary uncertainty is unavailable
where appropriate. Failed CV folds and nonconverged power-study fits remain
visible instead of being silently treated as successful estimates.

## Remaining release gates

1. Run the plan's literal interactive RStudio check, or explicitly accept the
   command-line equivalent. The command-line checks
   rebuild vignettes using RStudio's Pandoc, but no GUI-session check is claimed.
   Source validate_in_rstudio.R from the project root to save that receipt.
2. Publish the required landgraph version before submitting terradish to CRAN.
3. Obtain owner sign-off under plan Section 12 before setting version 0.1.0,
   merging to master, tagging, or running the release synchronization.

The owner's downstream manuscript and SLiM refits are expressly outside this
package task. RESUME.md identifies the next steps and saved validation artifacts.

Final diff cleanup corrected two obsolete setup comments in the spline and
model-comparison vignette sources. This comment-only change followed the final
archive check; evaluation settings, executable chunks, and rendered prose are
unchanged. The archive fingerprint refers to the preserved checked candidate.
