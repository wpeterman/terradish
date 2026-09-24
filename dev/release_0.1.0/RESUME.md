# terradish 0.1.0: resume checkpoint

Updated September 24, 2026 after final local validation.
Implementation is committed at development version 0.0.57 on release/0.1.0.
Use git log -1 for the final commit hash. RELEASE_REVIEW.md is the owner-facing
review; log.md contains phase and numerical evidence. The original logoff
checkpoint is preserved in logoff_status.md and Git history.

## Completed work

All eight implementation/validation phases are complete locally. The final
candidate is dev/check/phase8-final-package. It has been adopted into this
checkout. Do not adopt earlier candidates or copy compiled artifacts.

- Clean source tests: 988 passing assertions, no failures/errors, 197 diagnostic
  warnings, one intentional parallel skip.
- Full .57 archive check: 0 errors, 0 warnings, 0 notes; 990 installed assertions
  pass with 197 diagnostic warnings and no skips. Evidence: phase8_check.*.
- Final guide/caption-corrected archive check: 0 errors, 0 warnings, 0 notes.
  Examples and all eight guides rebuilt. Tests were skipped because R, C++, and
  test sources match the complete passing .57 check byte for byte.
  Evidence: phase8_final_check.*; archive: dev/check/phase8-final/terradish_0.0.57.tar.gz.
- Final archive SHA-256:
  dade1f8982f39ae55cf9e5f8a5b32cc87e0c3d55c5140cb74801c4797981d5cc.
  See phase8_final_fingerprint.json and phase8_final_vignette_archive.log.
- R line coverage: 86.79% overall, 88.03% across Phases 2-5 changed files.
  Three individual files remain below 80%. Native gcov collection failed with
  exit code 6; no native line-coverage claim. See phase8_r_coverage.* and
  phase8_validation_environment.md for successful and interrupted run boundaries.
- All nine numerical baselines converge and meet corrected C2 tolerances.
  Audits and their adaptations are preserved in phase8_audits and
  phase8_audit_interpretation.md. Fixed CV selects the noise extension in 3/8
  Wishart and 2/8 MLPE replicates, excluding one numerical tie.
- Fast CRAN numerical subset: eight assertions passed in 1.21 seconds combined.

## Branches and prerequisites

- release/0.1.0: final core implementation, version 0.0.57.
- experimental: 38f57d8, version 0.0.48; implementation unchanged from original
  0.0.47, with known limitations documented. Worktree is the sibling
  ../terradish-experimental-release-notes.
- landgraph: codex/terradish-core-companion, 36a4e17, version 0.0.3; local
  checks pass. Private validation library contains this version.
- Original master 6a2c13f and v0.0.47 remain intact. No push, merge to master,
  new release tag, or CRAN submission has been performed.

## Next actions

1. Review RELEASE_REVIEW.md. The plan requests a literal interactive RStudio
   check. Run source("dev/release_0.1.0/validate_in_rstudio.R") in RStudio from
   the project root, or obtain explicit owner acceptance of the completed
   command-line equivalent. No GUI-session execution is claimed.
2. Obtain owner sign-off before version 0.1.0, master merge, tag, push, or release
   synchronization, as required by plan.md Section 12. All requested local
   commits were already authorized; do not ask again for that authority.
3. Publish the required landgraph version before terradish CRAN submission.
   Do not treat local checks as CRAN acceptance or publication evidence.
4. If releasing, update versioned metadata and rerun validation for the actual
   release artifact. Preserve existing evidence rather than overwriting it.

## Owner decisions already authorized

Preserve the original vignette rewrapping; use truthful covariance-centering
metadata; make Gower the grouped-data default; investigate/fix Gaussian
stopping numerics without relaxing the specified thresholds. All implemented.
Downstream manuscript and SLiM refits remain outside this package task.

Final diff cleanup corrected two obsolete setup comments in the spline and
model-comparison vignette sources. This comment-only change followed the final
archive check; evaluation settings, executable chunks, and rendered prose are
unchanged. The archive fingerprint refers to the preserved checked candidate.
