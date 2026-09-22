## Submission summary

This is a new submission of terradish, a maximum-likelihood estimator for
parameterized conductance (resistance) surfaces in landscape genetics.

## Test environments

* Local Windows 11, R 4.6.1, checked from the exact source archive:
  0 errors, 0 warnings, 1 note.
* R-hub Ubuntu release: Status OK.
* R-hub Windows R-devel: examples, tests, and vignette rebuilding passed. The
  only warning was a CRLF line ending in `configure.ac` introduced by R-hub's
  Git checkout. The submitted source archive was independently verified to
  contain LF line endings and passed the local Windows check without warnings.
* win-builder R-release and R-devel: submitted on 2026-09-10; emailed results
  are pending.

## R CMD check results

The exact source archive produced 0 errors, 0 warnings, and 1 note under
`R CMD check --as-cran` on the local Windows environment. The one expected note
is the "New submission" note.

## Dependencies

* `landgraph` (>= 0.0.1) is available from CRAN; terradish's genetic-summary
  helpers re-export it.
* `multiScaleR` (Imports) is on CRAN.
* All other Imports are CRAN or base packages; Suggests are used conditionally.

## Additional checks

* `spelling::spell_check_package()` reported no spelling errors.
* The PDF and HTML manuals built successfully from the exact source archive.
* `urlchecker::url_check()` found only transient HTTP 504 responses from Zenodo
  for DOI 10.5281/zenodo.21225712. The DataCite registry reports that DOI as
  findable and registered to Zenodo.
