## Submission status

This is a corrected first submission of terradish 1.0.0. The initial archive
was submitted on 2026-10-05. CRAN incoming pre-tests on Windows and Debian
reported one NOTE: "New submission" and two possibly misspelled words in
DESCRIPTION, "multigrid" and "natively". Both words have been removed from
the DESCRIPTION prose without changing package code or numerical behavior.
The version remains 1.0.0 because the first submission was not accepted.

The required companion landgraph 0.0.3 is on CRAN, as is multiScaleR. The
1.0.0 release focuses on symmetric resistance models; research prototypes
remain on the experimental branch.

## Test environments

- Windows 11 x86_64, R 4.6.1 (R-release): full R CMD check --as-cran on the
  corrected submission archive, 0 errors, 0 warnings, 1 NOTE
  (New submission only).
- Windows 11 x86_64, R 4.5.3 (R-oldrelease): full R CMD check --as-cran on
  an earlier release metadata candidate, 0 errors, 0 warnings, 1 NOTE
  (New submission). This older check does not cover the corrected archive.

The corrected archive is terradish_1.0.0.tar.gz (SHA-256
6a0504037acd89e801f2f0de477916f61492ae4f1dd45293864d58709d56bc58,
2,974,270 bytes). The exact archive passed installation, ordinary and
donttest examples, installed tests, all vignette rebuilds, and PDF and
HTML manual checks. Installed tests reported 814 passing expectations, no
failures, 190 captured diagnostic warnings, and five intentional CRAN skips.
The diagnostic warnings did not produce package-check warnings.

The R 4.5.3 check used an earlier archive (SHA-256
224645503742cf2ab9811ccfad0186f308b0642b8b0fb7184d524134016adca9),
before the release date and Zenodo citation were corrected. GitHub Actions
passed all six Windows, macOS, and Ubuntu jobs on the released commit
fe2e4b7. Those CI results predate this DESCRIPTION-only correction. Local
checks are not CRAN acceptance.

## Note explanation

The only NOTE on the corrected archive is expected for a first submission:
"New submission". The two spelling flags from the initial incoming pre-tests
no longer occur in the corrected archive's R 4.6.1 incoming-feasibility check.

## Reverse dependencies

This is a first submission, so there are no CRAN reverse dependencies.
