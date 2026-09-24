## Submission status

This is preparation for the first CRAN release of terradish 0.1.0.
The current development version is 0.0.56. No submission has been made for
this core overhaul; owner release approval remains pending.

## Validation

Windows, R 4.6.1: the .54 robustness candidate passed a complete check with
0 errors, 0 warnings, and 0 notes. The .55 function-documentation candidate
passed examples and vignette rebuilding with the same clean result; its code
was unchanged from .54 and its redundant test run was skipped.
Final .56 source-archive and as-cran checks are pending. Earlier R-hub and
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
