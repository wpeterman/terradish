# Logoff status: 2026-09-24 14:42 America/New_York

Both remaining validation jobs were deliberately stopped before the owner's
14:45 logoff. These are interrupted checks, not successful full checks.

| Check | Completed | Interrupted at |
|---|---|---|
| C3 frozen .52 candidate | Installation, static checks, ordinary and extended examples, 897 installed-test assertions, vignette dependency checks | Rebuilding vignette outputs |
| C4 frozen .53 candidate | Installation, static checks, ordinary and extended examples | Installed-package tests |

Both frozen candidates reported the same NOTE about unqualified qchisq/tail.
The exact current source already fixes this using stats::qchisq and utils::tail.
That correction passed 11 focused profile assertions and an adopted-package
unresolved-global check. The full exact-current check is still required.

C3 installed tests: 897 passed, 0 failures, 40 warnings, no skips.
C4 broad source tests: 926 passed, 0 failures/errors, 194 warnings, one skip.
Neither warning count is represented as zero in the validation record.

Final logs:
- phase3_interrupted_00check.log
- phase3_installed_test.log
- phase4_interrupted_00check.log
- phase4_interrupted_installed_test.log (partial only)
- phase3_accepted_check.log and phase4_check.log retain the complete console
  output up to interruption. No completed-check RDS receipt was produced.

C3 process tree stopped at 14:41:52; C4 stopped at 14:41:11. No scheduled or
background continuation was created. All source changes, small validation
receipts, scripts, logs, and resume instructions are saved in the local C4
commit. Three large raw diagnostic fit objects remain local and Git-ignored.

Next action: follow RESUME.md and run validate_current_checkout.R with a fresh
label to complete the exact .53 vignette-enabled check before Phase 5.
