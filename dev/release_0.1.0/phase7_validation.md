# C7 documentation and examples

The .56 candidate rewrites README and all eight retained vignettes around the
core API. Its 45 R source files have identical executable syntax to C6. No
new numerical behavior is introduced in this documentation batch.

All eight vignettes rendered from source. Six receipts are in
`phase7_render_first.rds`; the IBE and simulation guides are in
`phase7_render_remaining_v2.rds`. The first IBE render exposed an example-only
integer/numeric `vapply` mismatch, fixed before the second render. Rendered
text and computational outputs were inspected; no claim of browser-based
visual inspection is made. Later executable clamping examples are covered by
the final package check.

The README's base and joint fits converge. The new IBE example's four models
converge under both likelihood families. Roxygen and package/vignette spelling
checks pass. Prose was checked with the suite's humanizer audit; its remaining
flags for "leverage" refer to the statistical term and are intentional.

The optional precompute script is retained as an explicit long-running recipe.
No active vignette reads its old saved outputs, and no old output is presented
as newly regenerated evidence. The radish benchmark script was inspected and
updated for convergence reporting and core settings; its timings were not
rerun and are not used as performance claims.

## Dependency and metadata review

The reviewed Imports remain used: ggplot2 for plots; MASS for a generalized
inverse fallback; nlme for least-squares GLS; multiScaleR for the external
kernel-scale workflow; parallel for worker management; splines for fitted
bases; utils for object sizes and tabular helpers; grDevices for colors; and
methods for class construction/coercion. No new dependency was added.

DESCRIPTION uses the core scope and requires landgraph >= 0.0.3. NEWS records
the breaking changes and the planned release. CITATION.cff, .zenodo.json,
cran-comments.md, and WORDLIST were updated. The existing title and DOI are
preserved. The final 0.1.0 version and release actions await owner sign-off.

## Scope scan

The literal Phase 1b grep has matches because `radish_cv` occurs inside retained
`terradish_cv_*` names. It also finds the shared compiled block-CG helpers
explicitly permitted by S9 and a test of the rejected `exact_refine = FALSE`
option. `phase7_scope_classified.log` records every exception. No public
NAMESPACE export, dispatch, usage, or documentation for moved features remains.

## Final validation

The .56 examples/vignette check passed with zero errors, warnings, or notes in
16m44s (`phase7_check.rds`). All three standalone example functions completed.
The spline examples' fits converge. The power example reports fit rates from
0.8 to 1, conductance-recovery power from 0 to 0.2, and an empty planning table:
none of these small designs reaches its 80% recovery target. This is an
illustration of the workflow, not a validated sample-size recommendation.
Its default run took about half an hour on this concurrently loaded machine.
Tests are
intentionally skipped in the .56 check because its executable R code is
identical to C5/C6; C8 runs the complete source and installed test suites.
