# terradish 0.1.0 implementation checkpoint

## Authorization

The owner authorized the supplied plan and local commits. The owner approved
including the original model-comparison paragraph wrapping in C0 and truthful
landgraph covariance metadata instead of unconditional centered="sites".
Do not push or merge master before owner sign-off. Stop active jobs and save
this checkpoint before 14:45 America/New_York on 2026-09-24.

## State at 12:40

- Original: 6a2c13f, version 0.0.47; preserved by v0.0.47 and experimental.
- Active branch: release/0.1.0. No terradish phase commits yet.
- IMPORTANT: Git index contains ONLY the accepted-by-tests Phase 1a helper
  relocation. Working tree contains the separate Phase 1b feature removals.
  Do not stage all files until C0 and C1a have been committed separately.
- C0 baseline: all nine models saved; source tests 920 pass, 0 failures/errors,
  one expected warning, one intentional skip. Full check passed installed
  tests and is rebuilding vignettes. Await final check result before C0.
- C1a: 414 function bodies/formals/classes unchanged; no new unresolved globals;
  all nine numerical baselines identical; source tests same as original.
  Full package check is running from its frozen temporary source copy.
- C1b: removals and generated docs prepared. Source tests, full vignette-enabled
  check, and nine-model rerun are running. Globals scan clean. Acceptance pending.
- Phases 2-9 not implemented. Mathematical Wishart contrast mapping and original
  Gaussian alignment reproduced in scripts/logs, ready for Phase 2.
- landgraph companion: L1/L2 committed as 501aa23, version 0.0.2, branch
  codex/terradish-core-companion. 73 tests and full check passed. L3 pending.
  Nothing installed, pushed or merged. terradish uses installed landgraph 0.0.1.

## Commit separation

C0 must contain preparation records and the original owner vignette only,
with version 0.0.48. The exact owner vignette is backed up at
`dev/check/owner-model-comparison.Rmd`. Back up the CURRENT split vignette
before temporarily restoring that owner version for C0. Use explicit paths
with git commit --only to preserve the staged C1a refactor. Restore the split
vignette afterward. Do not include unvalidated feature removals in C0.

Then bump to 0.0.49, stage only metadata and C1a records, and commit the
already-staged C1a source snapshots after its full check passes. Do not
restage R or man paths from the Phase 1b working tree into this commit.
C1b version will be 0.0.50 after its own validation passes.

## Validation

Use R 4.6.1 at C:/Program Files/R/R-4.6.1/bin/Rscript.exe.
Clear LC_* environment values equal to C.UTF-8. Set RSTUDIO_PANDOC to
C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools.
The validate_phase.R runner stages source without .git to avoid Windows
path-length failures. Do not edit scripts while they are running.

Each phase has *_test.rds, *_check.rds and logs in this directory. Full check
logs are also under dev/check/<phase>/terradish.Rcheck/00check.log.
Baseline originals must never be overwritten. core_baseline.R takes a new
output filename; compare_baseline.R compares it to the frozen original.

The plan's literal grep has two unavoidable false positives: radish_cv is a
substring of retained terradish_cv_folds; S9 explicitly permits internal
compiled solvers sharing radish.cpp. Negative tests also mention rejected
exact_refine = FALSE. Check exported names and call paths, not literal hits alone.

Original audit evidence: sibling Research/R_packages/terradish checkout,
review_work_20260924/package_audit. plan.md is the complete owner specification.
Use git --no-optional-locks for read-only Git commands.
