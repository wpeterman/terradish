# terradish 0.1.0: resume checkpoint

Saved 2026-09-24 before the owner's 14:45 America/New_York logoff.
Read logoff_status.md for the exact final check states and stopped jobs.
The overhaul is NOT complete. Phases 5-9 and final release validation remain.

## Current checkout and commits

Branch: release/0.1.0. Current development version: 0.0.53 (C4).
The C4 commit includes this checkpoint; use git log -1 for its hash.

- C0 c80425d, .48: baseline evidence and owner vignette wrapping.
- C1a 488ac6d, .49: helper relocation; 414 function contracts unchanged.
- C1b 604a3a1, .50: supported-core split (amended authoritative hash).
- C2 c3ac708, .51: F1-F8 numerical/default corrections; full check clean.
- C3 c2aa209, .52: inference and constrained-optimization repairs.
- C4 current HEAD, .53: per-model/repeated CV, fixed nuisance prediction,
  checkpoint signatures, failure accounting, and paired comparisons.

Original master 6a2c13f (.47) is preserved by v0.0.47 and experimental.
Nothing has been pushed, merged to master, or submitted for release.
Owner approval is still required for those release actions.

## Owner decisions already authorized

The full plan and local commits are authorized. Do not ask again about:
- Including the original vignette paragraph rewrapping in preparation.
- Truthful landgraph covariance-centering metadata.
- Gower as the grouped-data default.
- Investigating/fixing the stalled Gaussian baseline without relaxing the
  specified convergence thresholds. Those thresholds remain unchanged.

Exact plan: plan.md here; original Downloads/terradish_core_release_plan_0.1.0.md.

## Validation completed

C0/C1a/C1b/C2 full vignette-enabled checks passed. C2 had zero errors, warnings,
and notes; 803 source assertions passed. All nine C2 contracts were verified.

C3:
- 893 broad source assertions passed, zero failures/errors, 40 warnings and one
  intentional parallel skip (phase3_repaired_test.*).
- 44 final inference/guard assertions passed, zero failures/errors; 21 inner
  nuisance iteration-limit warnings (phase3_final_guard_test.*).
- All nine baseline fits converge (code 0), matching C2 estimates, nuisance
  parameters, SEs and log likelihoods within specified tolerances.
- Installed-package tests: 897 pass, 0 fail, 40 warnings, 0 skips
  (phase3_installed_test.log). Ordinary and extended examples passed.
- The full frozen check reported one NOTE for unqualified qchisq/tail.
  These are fixed in current source as stats::qchisq and utils::tail.
  The qualified helper passes 11 focused assertions without warnings and the
  adopted package has no added unresolved globals. No full clean-check claim
  is made until the exact current checkout completes its check.

C4:
- 926 broad source assertions passed, zero failures/errors, 194 warnings and
  one intentional parallel skip (phase4_test.*).
- 48 focused CV assertions passed, zero failures/errors (154 warnings, chiefly
  unused covariates plus recorded baseline/optimizer warnings).
- Noise covariate audit: Wishart extension preferred 3/8, MLPE 2/8.
- Real effect (lambda=1) audit: Wishart preferred 8/8, MLPE 5/8.
- Each audit excludes one stalled MLPE fold; that replicate retains three
  common folds. A failed model fold makes its full total NA.
- Nu=500 scores equal 20 times nu=25 scores with the same ranking.
- Documentation regenerated; globals and package-source diff checks pass.
- Current R/man/tests equal the frozen .53 candidate except the separately
  validated profile-helper qualifications (candidate_differences.txt records).

Warnings are retained in the logs; passing assertions are not a claim that
all simulated fits converge. See phase3_phase4_validation.md for interpretation.

## Numerical and CV repairs worth preserving

Newton steps now invert only the free Hessian block at active bounds.
Extreme MLPE rho warm starts are reset without changing the other warm values;
fixed predictive parameters are exempt. Stable correlation subspace algebra
avoids upper-bound Woodbury cancellation. Joint nuisance covariance includes
conductance uncertainty. A zero lambda is separate from no spatial structure.

Fixed-nuisance CV can score no-structure fits because beta/tau=0 makes their
predictions independent of conductance. Reprofiled unidentified fits still fail.
Failed baseline scores are separate from valid trained-model scores; their
corresponding gains stay NA. Rankings use folds where every model succeeded.

## First actions on resume

1. Read logoff_status.md and git status. Do not rerun already-passed source tests
   without a code change or new concern. Complete the full vignette-enabled
   check on the exact current checkout, including the namespace fix:

   Rscript dev/release_0.1.0/validate_current_checkout.R phase4_exact

   Run from the package root with the Windows environment setup below. Redirect
   output to dev/release_0.1.0/phase4_exact.log. The runner refuses to overwrite
   existing evidence; choose a new label for a later rerun.
2. Resolve any failures, then implement Phase 5 (C5, .54), Phase 6 (.55),
   Phase 7 (.56), and remaining .1.0 release/experimental-branch gates.
3. During Phase 6, replace the CV help example's maxit=1 with a converged,
   lightweight example; it now correctly demonstrates failed training fits.
4. Do not push, merge master, or publish without owner sign-off.

## Companion, environment, and saved artifacts

../landgraph branch codex/terradish-core-companion is clean:
501aa23 (.2 metadata), 36a4e17 (.3 gower default). Both passed full checks.
Gower recovered expected signs in all three L3 seeds; within diagonal did not.
Companion publication remains an external dependency gate.

Use C:/Program Files/R/R-4.6.1/bin/Rscript.exe. Before launching R, clear LC_*
environment values equal to C.UTF-8. Set RSTUDIO_PANDOC to
C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools.
Use .libPaths(c("dev/check/phase2-library", .libPaths())): that private library
contains landgraph .3; the ordinary installed library still has .1.

Frozen candidates and check directories remain under dev/check; their evidence
must not be overwritten. Editable prototype snapshots are also saved in
work_in_progress/phase3 and phase4 here. Main source now contains both phases.
The three large raw diagnostic RDS files are retained locally and explicitly
ignored by Git; their scripts/logs and acceptance receipts are tracked.
Read-only Git uses --no-optional-locks. Raw logs and source snapshots preserve
original whitespace. Never delete original baselines or unrelated scratch work.
