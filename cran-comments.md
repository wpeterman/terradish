## Submission status

This is preparation for the first CRAN release of terradish 0.1.0.
The current development version is 0.0.60. No submission has been made for
this core overhaul; owner release approval remains pending.

## Validation

Version 0.0.60 binds Gaussian factory settings at construction so factories
built in a loop retain their own raster choices and scale settings. The new
14-expectation regression failed seven expectations on 0.0.59 and passed on
0.0.60. The full source test suite completed without failures or errors, with
198 diagnostic warnings and one intentional parallel skip. The installed
archive tests passed 806 assertions with no failures, 190 diagnostic warnings,
and five expected skips. Objective, gradient, Hessian, exact and coarse-raster
fits, and prediction were identical between the two builds for a concrete
Gaussian factory call. The vignette-enabled source archive
`terradish_0.0.60.tar.gz` (SHA-256
`218d9a86114740c08beed4dbfc8ff1eef21aa8cce5fdc6eff0976051af151c23`)
passed `R CMD check --no-manual` on Windows R 4.6.1 with Status: OK.

Version 0.0.59 fixes Gaussian map-unit widths in
`covariance_response_power()`. The focused Gaussian refit regression passed.
`devtools::document()` regenerated the affected help page, and the full source
test suite passed with 991 assertions, 0 failures, 197 expected diagnostic
warnings, and one intentional parallel skip. The vignette-enabled source archive
`terradish_0.0.59.tar.gz` (SHA-256
`6b8597a1c1bae7e8a99537c1aeeed12a6a0ba50dfebf7f9ff9b4c707312a5bdb`)
passed `R CMD check --no-manual` on Windows R 4.6.1 with Status: OK, including
the installed tests (792 passing assertions, 190 diagnostic warnings, five
expected skips), examples, and all eight vignette rebuilds. The corrected W1
run later completed with 720 fitted cells, 80 structural skips, and no top-level
cell errors. Its scientific interpretation remains under review.

The 0.0.58 follow-up changes documentation and metadata only. All eight
vignettes rebuilt successfully, all 62 Rd pages passed checkRd, and browser
inspection found no rendering errors in 61 vignette math expressions and 124
help-page math expressions. Executable R syntax is unchanged from 0.0.57.
The complete numerical test suite and package archive check were not repeated
for this documentation-only follow-up; the earlier evidence below remains
specific to the .57 archives.

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
