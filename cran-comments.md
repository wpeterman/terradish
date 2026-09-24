## Submission status

This is preparation for the first CRAN release of terradish 0.1.0.
The current development version is 0.0.57. No submission has been made for
this core overhaul; owner release approval remains pending.

## Validation

Windows, R 4.6.1: the complete .57 development archive passed
devtools::check(args = "--as-cran") with 0 errors, 0 warnings, and 0 notes,
including 990 passing installed test assertions, examples, donttest examples,
and all eight rebuilt vignettes. The final guide/caption-corrected archive
also passed with 0 errors, 0 warnings, and 0 notes. Its tests were skipped
because R, C++, and test sources are byte-identical to the complete passing
check. devtools used its no-manual default and disabled remote incoming
checks; this is local validation, not CRAN acceptance.

A separate clean source test run passed 988 assertions with no failures/errors,
197 diagnostic warnings, and one intentional parallel skip. Installed tests
had the same warnings and no skips. R-source coverage is 86.79% overall and
88.03% across changed files; native coverage collection failed and remains
unverified.

The command-line vignette check used RStudio's bundled Pandoc. The release
plan's literal interactive RStudio check remains open. Earlier R-hub and
win-builder runs predate this overhaul and do not validate this candidate.

## Dependencies

terradish requires landgraph >= 0.0.3 for the documented covariance metadata
and grouped-data Gower default. This companion version has passed local checks
but its CRAN publication is still a release prerequisite. The private validation
library contains that version. Do not submit terradish before the required
landgraph release is available on CRAN.

## Scope

The core retains symmetric resistance likelihoods, supported conductance
factories, direct/auto/AMG solvers, and fixed-domain cross-validation.
Research prototypes are retained on the experimental branch. Final check
results and submission metadata will be updated after release approval.
