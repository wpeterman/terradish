## Submission summary

This is a new submission of terradish, a maximum-likelihood estimator for
parameterized conductance (resistance) surfaces in landscape genetics.

## Test environments

* Local: Windows 11, R 4.6.0.
* (Before submission) win-builder: R-release and R-devel via
  `devtools::check_win_devel()`.
* (Before submission) R-hub across the standard CRAN platforms.

## R CMD check results

Target: 0 errors | 0 warnings | 1 note.

The one expected note is the "New submission" note. Please update this file with
the final `R CMD check --as-cran` result from win-builder and R-hub before
submitting.

## Dependencies

* `landgraph` (Imports) has been submitted to CRAN and should be available before
  or with this submission; terradish's genetic-summary helpers re-export it.
* `multiScaleR` (Imports) is on CRAN.
* All other Imports are CRAN or base packages; Suggests are used conditionally.

## Notes for the maintainer before submitting

* Run `spelling::spell_check_package()` and reconcile `inst/WORDLIST`.
* Run `urlchecker::url_check()` on the README and vignettes.
* Confirm the PDF manual builds on a clean TeX toolchain (a local MiKTeX
  `xkeyval` issue produced a spurious PDF-manual error here).
