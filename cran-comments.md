## Submission status

This is a corrected first submission of terradish 1.0.0. The initial archive
was submitted on 2026-10-05. CRAN incoming checks flagged two possibly
misspelled words in DESCRIPTION, "multigrid" and "natively". Both have been
removed without changing package code or numerical behavior. The version
remains 1.0.0 because the first submission was not accepted.

The CRAN Windows check also took 19 minutes, with tests taking 392 seconds
and vignette rebuilding 400 seconds. Routine checks now use small, explicitly
labeled teaching subsets in the computational vignettes. Expensive parallel,
simulation, and full-resolution regression tests are conditional on
`NOT_CRAN=true`; the GitHub Actions workflow sets that variable so those tests
still run on Windows, macOS, and Linux.

The required companion landgraph 0.0.3 is on CRAN, as is multiScaleR. The
1.0.0 release focuses on symmetric resistance models; research prototypes
remain on the experimental branch.

## Test environments

- Windows 11 x86_64, R 4.6.1: full `R CMD check --as-cran` on the exact
  corrected source archive, 0 errors, 0 warnings, 1 NOTE (New submission).
  Installation, ordinary and donttest examples, installed tests, vignette
  rebuilding, and PDF and HTML manual checks passed.
- Windows 11 x86_64, R 4.5.3: an earlier release metadata candidate passed
  `R CMD check --as-cran` with the same expected NOTE. This does not cover
  the current archive.

The current archive is `terradish_1.0.0.tar.gz`, 2,641,784 bytes, SHA-256
`52d6dc542d56de7e87f17ddd7ef4452d2a323ef946eb2c64692c80d55b825605`.
In the same local R 4.6.1 setup, installed tests fell from 11 minutes to
201 seconds, and vignette rebuilding fell from 13 minutes to 230 seconds.
Installed tests reported 784 passing expectations, no failures, 187 captured
diagnostic warnings, and 10 intentional CRAN skips. The diagnostic warnings
did not produce package-check warnings. CRAN's own timing on this archive
has not yet been measured.

## Note explanation

The only NOTE on the corrected archive is expected for a first submission:
"New submission". The two spelling flags from the initial incoming checks
do not occur in the corrected archive's R 4.6.1 incoming-feasibility check.

## Reverse dependencies

This is a first submission, so there are no CRAN reverse dependencies.
