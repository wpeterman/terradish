## Development validation, October 28, 2026

This record describes the combined development tree on
`codex/amphibian-vignette`, not a submitted or accepted CRAN archive.
Version remains 1.0.0 pending an authorized commit and release decision.

## Incorporated review corrections

- Added the verified Clarke, Rothery and Raybould (2002)
  <doi:10.1198/108571102320> and Hanks and Hooten (2013)
  <doi:10.1080/01621459.2012.724647> method references to DESCRIPTION.
- Added return documentation for the nine legacy topics and the plot-list
  print method; regenerated help from the combined roxygen sources.
- Preserved worker session settings through explicit source paths and
  qualified external function calls. Source-loaded workers also load the
  new MLPE and Wishart profiling helpers.

## Additional development changes

The integration includes deterministic MLPE and Wishart nuisance profiles,
unit-aware nuisance curvature, direct-solve reference and residual checks,
fixed/custom-bound protections, and outer/final-nuisance convergence guards.
A failed nuisance warm start also retries the default values within the same
bounds, accepting only a converged retry. These numerical changes have
dedicated regression coverage. Public
warnings and convergence thresholds were not removed or relaxed.

The eleven-guide documentation suite includes pond examples, data-generation
provenance, categorical predictors, and the individual-level demonstration.
Conductance maps use consistent logarithmic viridis colors and labels in
conductance units. Estimate and interval maps share limits within a fit.
Computational guides retain precomputed displays for routine package checks;
their executable sources are distributed in inst/vignette-source.

## Validation of this combined tree

- Windows, R 4.6.1: roxygen documentation and installation completed.
- Extended installed-package suite: 216 tests, 1,803 passing expectations,
  zero failures, zero errors, and zero warnings. The one source-worker
  skip is covered by a separate source-loaded run.
- Source-loaded worker regressions: 14 passing expectations, zero failures,
  zero errors, and zero warnings. The installed-only case is covered by
  the extended installed suite.
- All seven computational-guide sources were executed against the combined
  package. All 23 saved fits and 42 cross-validation rows passed convergence
  and finiteness checks. Numerical differences from the earlier examples
  were negligible and do not change their interpretation.
- All eleven vignettes rendered. The HTML audit found no issues across 46
  embedded figures, their descriptions, local links, and MathML. All 16
  conductance panels were visually inspected after the final rebuild.
- Windows, R 4.6.1: a fresh vignette-enabled R CMD check --as-cran --no-manual
  completed with 0 errors, 0 warnings, and 1 NOTE (New submission). Ordinary
  and donttest examples, routine tests, and vignette rebuilding passed.
  The PDF manual was not built in this run.

No fresh macOS, Linux, R-devel, or CRAN-host timings are claimed for this
integration. Earlier archive results do not validate this development tree.
Detailed receipts and source hashes are retained under
dev/vignette_consistency_20261028.

Checked archive: terradish_1.0.0.tar.gz, 6410061 bytes.
SHA-256: 0d5a867bee8d0d0d175683822940fef0abe7bcc87825ba10f5129c7e6d1d6b38.

After this check, a development-history paragraph was removed from the
categorical vignette and that guide was rendered again. No code, data,
results, or figures changed; the package check predates this prose-only edit.
