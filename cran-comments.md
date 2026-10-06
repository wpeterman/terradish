## Submission status

This is a further corrected first submission of terradish 1.0.0. The initial
archive was submitted on 2026-10-05. CRAN incoming checks flagged two
possibly misspelled words in DESCRIPTION, "multigrid" and "natively". Both
were removed without changing package code or numerical behavior. The version
remains 1.0.0 because the first submission was not accepted.

The initial CRAN Windows check took 19 minutes, with tests taking 392 seconds
and vignette rebuilding 400 seconds. The next archive reduced these to 172
and 108 seconds on CRAN Windows, respectively, but its additional-issue check
still reported overall checktime of 11 minutes, above CRAN's 10-minute limit.
For this archive, more multi-fit integration and simulation regressions are
conditional on `NOT_CRAN=true`; the GitHub Actions workflow sets that variable
on Windows, macOS, and Linux. Routine CRAN checks retain 387 passing test
expectations. Computational vignettes continue to use small, explicitly
labeled teaching subsets. Package code and numerical behavior are unchanged.

The required companion landgraph 0.0.3 is on CRAN, as is multiScaleR. The
1.0.0 release focuses on symmetric resistance models; research prototypes
remain on the experimental branch.

## Test environments

- Windows 11 x86_64, R 4.6.1: full `R CMD check --as-cran` on this exact
  source archive, 0 errors, 0 warnings, 2 NOTEs (New submission and a local
  MiKTeX temporary file, `lastMiKTeXException`, left after the successful
  PDF manual check). Installation, ordinary and donttest examples, installed
  tests, vignette rebuilding, and PDF and HTML manual checks passed.
- Windows 11 x86_64, R 4.5.3: an earlier release metadata candidate passed
  `R CMD check --as-cran` with the same expected NOTE. This does not cover
  the current archive.

The current archive is `terradish_1.0.0.tar.gz`, 2,648,274 bytes, SHA-256
`a8bab3dd1413cd2b641f65dd16647b21da26c7d6c987d4dc9a05fc1bdf63462b`.
In the same local R 4.6.1 setup, installed tests fell from 201 to 57 seconds;
vignette rebuilding took 197 seconds. Installed tests reported 387 passing
expectations, no failures, one captured diagnostic warning, and 22 intentional
CRAN skips. The captured warning did not produce a package-check warning.
The full extended source suite passed with no failures before the archive was
built. CRAN's own checktime on this new archive has not yet been measured.

## Note explanation

"New submission" is expected for a first submission. The second local NOTE
is a MiKTeX temporary file created during PDF manual validation on this
Windows host, not by package code. The two DESCRIPTION spelling flags from
the initial incoming checks are absent.

## Reverse dependencies

This is a first submission, so there are no CRAN reverse dependencies.
