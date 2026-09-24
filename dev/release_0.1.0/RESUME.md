# terradish 0.1.0 implementation checkpoint

Updated 2026-09-24 13:20 America/New_York. Owner logs off at 14:45 today.
Save a final checkpoint and stop this task's active jobs before then.

## Authorization and decisions

Owner authorized the complete plan and local phase commits. No push or merge
master before owner sign-off. Owner approved original vignette paragraph
wrapping in C0, truthful covariance-centering metadata in landgraph, and
changing grouped-data genetic covariance to the gower default.

## Commits and separation

- Original 6a2c13f, version 0.0.47, preserved by v0.0.47 and experimental.
- Active branch release/0.1.0.
- C0 c80425d, version 0.0.48: preparation, baselines, owner vignette edit.
- C1a 488ac6d, version 0.0.49: behavior-preserving helper relocation.
- IMPORTANT: Git index now contains Phase 1b removals only. Worktree matches
  these staged package changes. Do not add Phase 2 before committing C1b.
- Phase 2 complete candidate is frozen in dev/check/phase2-package. Prototype
  inputs are in dev/check/phase2/R and tests; prepare_phase2.R records assembly.
- Phase 3 partial, unvalidated prototypes are in dev/check/phase3/R. Do not
  adopt them as a completed phase. I3 acceptance is on hold (see below).
- Phases 4-9 are not implemented. Full specification is plan.md.

## Validation

- Baseline: 920 assertions pass, zero failures/errors, one expected warning and
  one intentional skip. Full vignette-enabled check Status: OK (52m42s).
- C1a: all 414 function bodies/formals/classes unchanged; globals clean;
  nine numerical baselines identical. 920 assertions pass. Full check Status:
  OK, zero errors/warnings/notes (53m48s). phase1a_00check.log is durable.
- C1b: globals/scope checks pass; nine baselines identical; 703 assertions
  pass, zero failures/errors, one expected warning and one skip. Full check
  running, installed tests passed, rebuilding vignettes. Session 39852.
- Phase 2: 803 assertions pass, zero failures/errors, 11 warnings (10 new
  intentional Gaussian-support warnings and one legacy maxit warning), one
  skip. Result saved in phase2_test.rds. The runner has a trailing top-level
  else parse error AFTER saving its result. Fix the runner only after the
  concurrent check using that script finishes. Full check session 58503.
- Phase 2 prototype nine-model comparison passes the specified contracts;
  actual complete-package baseline is running, session 11996, core_phase2.rds.
- Phase 3 outer-optimizer prototype baseline: eight fits code 0; corrected
  even-grid Gaussian fit code 2 (stalled), projected gradient 0.0007226463.
  Owner informed and asked whether to investigate while retaining thresholds.
  Do not proceed past this explicit plan acceptance gate without resolution.
  This test binds the new outer optimizer but not the inner warm-start path;
  complete-package testing is still required for Phase 3.

## landgraph companion

Sibling ../landgraph, branch codex/terradish-core-companion:
- 501aa23, 0.0.2, metadata and documentation. 73 assertions and full check pass.
- 36a4e17, 0.0.3, owner-approved gower default. 74 assertions and full check pass.
- L3 three-seed comparison: within diagonal gives wrong x1 sign in all three;
  gower recovers both signs and keeps both coefficients within two SEs.
- Installed 0.0.3 only in dev/check/phase2-library. Ordinary library remains
  0.0.1. Candidate DESCRIPTION requires landgraph >= 0.0.3. External companion
  release remains a release gate. Nothing pushed or merged.

## Next steps

1. Finish C1b check; copy its 00check.log here; bump to 0.0.50 and update NEWS;
   commit staged split plus explicit phase records. Do not stage Phase 2 yet.
2. Finish Phase 2 check and actual-candidate baseline comparison. Adopt only
   candidate source/docs/tests/NAMESPACE and dependency change. Bump to 0.0.51,
   record all F1-F8 numerical changes in NEWS, and commit after acceptance.
3. Resolve Phase 3 stopping-rule gate. Complete I1-I12 and tests before C3.
   Current prototypes are unfinished, especially profile intervals, boundary
   inference, nesting checks, and IBE ratios. See code, not merely this list.
4. Continue Phases 4-9 in order. Final version is 0.1.0 only after release gates.

Use R 4.6.1: C:/Program Files/R/R-4.6.1/bin/Rscript.exe. Clear LC_* values
C.UTF-8; set RSTUDIO_PANDOC to
C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools. Use login:false.
Never edit a running script or frozen candidate. Checks stage source without
.git to avoid Windows long paths. Read-only Git uses --no-optional-locks.

All original baseline files must remain immutable. New numerical runs use a
new filename. Raw log whitespace is intentional; exclude *.log from diff
whitespace checks. Plan literal grep has permitted false positives in internal
compiled solvers, retained terradish_cv_folds, and rejection tests.
