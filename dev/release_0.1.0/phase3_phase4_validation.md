# Phase 3 and 4 validation record

This record distinguishes completed tests from package checks still running.
The main package now contains C3 and C4 as development version 0.0.53. Final package-check states are recorded in logoff_status.md.

## Phase 3 inference and numerical repairs

All nine complete-package baseline fits converge with code 0. Their estimates,
nuisance parameters, standard errors, and log likelihoods match Phase 2 within
the plan's tolerances. The new stopping thresholds were not relaxed.

The numerical causes were an incorrect Hessian inverse at active bounds and
extreme MLPE correlation warm starts that could stay at a false zero-correlation
optimum. A separate predictive audit exposed subtractive cancellation near the
upper correlation limit. The repaired implementation uses the free Hessian
block and stable correlation subspace algebra.

Completed checks:

- Broad source suite: 893 assertions passed, zero failures/errors, 40 warnings,
  one intentional parallel skip (phase3_repaired_test.rds).
- Final nesting guard: 44 assertions passed, zero failures/errors, 21 inner
  nuisance iteration-limit warnings (phase3_final_guard_test.rds).
- Qualified profile helper: 11 assertions passed, no warnings/failures/errors;
  no undefined function references (phase3_namespace_fix_test.rds).
- Dedicated correlation tests compare dense inverses, determinants, quadratic
  forms, and extreme-correlation likelihoods (test-mlpe-boundary-stability.R).
- The joint nuisance covariance and IBE-ratio delta-method uncertainty agree
  with numerical joint-Hessian inversion (test-core-inference.R).

The full .52 candidate check is recorded in phase3_accepted_check.*. It was
frozen before stats::qchisq and utils::tail were qualified; its known NOTE must
not be reported as a clean check. The qualified helper is in work_in_progress.
Read the completed receipt, or RESUME.md if interrupted, before claiming a pass.

## Phase 4 cross-validation

Final focused contracts: 48 assertions passed, zero failures/errors, 154
warnings. Most warnings concern intentionally unused raster covariates; the
log also records convergence and baseline failures. Warnings were retained.
Tests verify fixed training nuisance parameters, mixed conductance factories,
repeated folds, checkpoint rejection after a changed response, paired-score
uncertainty, valid no-structure predictions, and failed-fold handling.

The nu=500 scores equal 20 times nu=25 scores with the same rankings. A failed
uniform-baseline optimization is recorded separately and cannot erase a valid
trained-model score; its gain remains unavailable.

Predictive audits use the supplied eight noise-audit seeds, random four-fold
partitions, a common landscape graph, and fixed training nuisance parameters.
Scores are summed only over folds where both models succeed. Full totals are
NA whenever a model fold fails. Numerical ties are zero only within 100 machine
epsilons of score scale.

| Audit | Wishart prefers extension | MLPE prefers extension | Failed folds |
|---|---:|---:|---|
| Noise covariate | 3 of 8 | 2 of 8 | One stalled MLPE fold; three common folds in that replicate |
| Environmental effect, lambda=1 | 8 of 8 | 5 of 8 | One stalled MLPE fold; three common folds in that replicate |

Noise results: phase4_predictive_final.rds/.log, frozen Phase 4 candidate.
Effect results: phase4_predictive_boundary.rds/.log, same predictive logic;
later edits only separate baseline errors, and these audits use baseline=FALSE.
The effect simulation uses nu=25 and S ~ Wishart(E + K + 0.1 I)/25.
These are acceptance simulations, not a general power or calibration claim.

Earlier effect runs failed because valid no-structure fits had no reported
conductance coefficients and were rejected. With beta/tau fixed at zero,
predictions are independent of those coefficients. A factory default now
supplies them for fixed predictive scoring only. The regression test confirms
that varying those coefficients leaves the score unchanged.

Full Phase 4 source and vignette-enabled package checks are recorded in
phase4_test.* and phase4_check.*. Consult RESUME.md for final completion or
interruption status. The frozen candidate shares the same already-fixed
profile-helper namespace NOTE described above.
