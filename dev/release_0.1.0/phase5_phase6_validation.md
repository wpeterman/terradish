# Phase 5 and Phase 6 validation

2026-09-24. Candidates were frozen during their checks and checked without
repository metadata to avoid Windows deep-path copying warnings.

## Acceptance evidence

- Exact C4 archive: complete check, 0 errors, 0 warnings, 0 notes.
  `phase4_exact_archive_check.rds` and matching log.
- C5 (.54): complete check, including examples, run-donttest examples,
  installed tests, and rebuilt vignettes. 0 errors, 0 warnings, 0 notes.
  Duration 21m12s; installed tests 562s; vignette rebuild 333s.
  `phase5_check.rds` and matching log.
- All nine C5 baseline fits converge (code 0) and meet the plan's tolerances
  against C2. Three duplicate-cell diagnostics are intentional new notices.
  `core_phase5.rds`, `phase5_baseline.R`, and matching log.
- Final affected C5 source tests: 167 passing assertions, no failures or
  errors, 154 diagnostic warnings. The complete installed test run above is
  the final whole-package acceptance evidence.
- C6 (.55): 45 R files have identical executable expressions to C5, verified
  by `phase6_source_contract.R`. Roxygen and package spelling pass. The
  documentation check ran examples, run-donttest examples, and vignette
  rebuilds; tests were deliberately skipped because C5 checks the same code.
  0 errors, 0 warnings, 0 notes; duration 12m27.9s.
- Prose audit and manual editing removed flagged wording from the new function
  documentation. The vignette/README rewrite belongs to C7 and remains separate.

## Development failures and their resolution

The initial broad C5 test run had 971 passing assertions, no assertion failures,
and three errors: two CV checkpoint resumes and one legacy factory-error text
expectation. Stored prediction metadata had captured a temporary formula
environment. Signature construction now strips environment attributes in that
metadata while retaining the fitted design and recipe. The legacy error keeps
its diagnostic phrase. Focused reruns and the complete check pass.

An initial robustness fixture accessed `$S` on a simulation result, which R
partially matched to `$Sigma`, the noiseless mean. That produced an exact MLPE
fit without a finite residual precision optimum. The fixture now uses the
sampled `$covariance`. An attempted Newton adjustment was fully reverted;
the C4 stopping thresholds and Newton implementation were preserved.

The separate exploratory `phase5_power_audit` predates the final convergence
guard. Do not report its counts as final acceptance. The final packaged
30-replicate nu_fit regression test passes with the guard applied.

Interrupted diagnostic and check logs are retained with their original names
or explicit interrupted suffixes. They do not supersede the final receipts.

## Remaining release work

C7 vignette and metadata revision, Phase 8 final audits/coverage/as-cran check,
and Phase 9 experimental-branch documentation remain. No release, master merge,
push, or companion publication has been performed.
