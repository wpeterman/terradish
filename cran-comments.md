## Submission status

This is a further corrected first submission of terradish 1.0.0. The initial
archive was submitted on 2026-10-05. CRAN incoming checks flagged two
possibly misspelled words in DESCRIPTION, "multigrid" and "natively". Both
were removed without changing package code or numerical behavior. The version
remains 1.0.0 because the first submission was not accepted.

The initial CRAN Windows check took 19 minutes, with tests taking 392 seconds
and vignette rebuilding 400 seconds. The next two archives reduced those
phases, but CRAN still reported overall checktime of 11 and then 15 minutes.
In the latest check, tests took 76 seconds and vignette rebuilding took 166
seconds; the remaining check stages and installation accounted for most of
the total. This archive retains the full prose, code, printed results, and
figures in five computational guides as precomputed displays. Their fully
executable R Markdown sources are included under `inst/vignette-source/` and
are rebuilt in continuous integration. Three other guides continue to execute
during routine vignette checks. The longer multi-fit integration and
simulation regressions remain conditional on `NOT_CRAN=true` and run in the
GitHub Actions matrix. Package code and numerical behavior are unchanged.

The required companion landgraph 0.0.3 is on CRAN, as is multiScaleR. The
1.0.0 release focuses on symmetric resistance models; research prototypes
remain on the experimental branch.

## Test environments

- Windows 11 x86_64, R 4.6.1: full `R CMD check --as-cran` on this exact
  source archive, 0 errors, 0 warnings, 1 NOTE (New submission). Installation,
  ordinary and donttest examples, installed tests, vignette rebuilding, and
  PDF and HTML manual checks passed.
- Windows 11 x86_64, R 4.5.3: an earlier release metadata candidate passed
  `R CMD check --as-cran` with the same expected NOTE. This does not cover
  the current archive.

The current archive is `terradish_1.0.0.tar.gz`, 2,948,822 bytes, SHA-256
`033fb5061663a8c1544c7ad531bcd330e904226028363a619654615e9413dd92`.
In the local R 4.6.1 check, installed tests took 29 seconds and vignette
rebuilding took 38 seconds. In the previous local archive, these phases took
57 and 197 seconds. Installed tests reported 387 passing
expectations, no failures, one captured diagnostic warning, and 22 intentional
CRAN skips. The captured warning did not produce a package-check warning.
The full executable source of each precomputed guide was rebuilt to create
the displayed results. CRAN's own checktime on this new archive has not yet
been measured.

## Note explanation

"New submission" is expected for a first submission. The two DESCRIPTION
spelling flags from the initial incoming checks are absent.

## Reverse dependencies

This is a first submission, so there are no CRAN reverse dependencies.
