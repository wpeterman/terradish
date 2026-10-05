## Submission status

This is the first CRAN submission candidate for terradish 1.0.0. It has not
been uploaded to CRAN. The required companion landgraph 0.0.3 is on CRAN, as
is multiScaleR. The 1.0.0 release focuses on symmetric resistance models;
research prototypes remain on the experimental branch.

## Test environments

- Windows 11 x86_64, R 4.6.1 (R-release): full `R CMD check --as-cran`
  on the submission archive, 0 errors, 0 warnings, 1 NOTE (New submission).
- Windows 11 x86_64, R 4.5.3 (R-oldrelease): full `R CMD check --as-cran`
  on the preceding metadata candidate, 0 errors, 0 warnings, 1 NOTE
  (New submission).

The submission archive is `terradish_1.0.0.tar.gz` (SHA-256
`0d12f2b8b866a21a157f05f31fe09f891235e29701bcd3848e01123b34b011be`,
2,975,048 bytes). The R 4.6.1 check ran ordinary and `\donttest{}` examples,
installed tests, all vignette rebuilds, and the PDF and HTML manuals.
Installed tests reported 814 passing expectations, no failures, 190 captured
diagnostic warnings, and five intentional CRAN skips. The warnings did not
produce package-check warnings. The R 4.5.3 check used the earlier archive
(SHA-256 `224645503742cf2ab9811ccfad0186f308b0642b8b0fb7184d524134016adca9`),
before the release date and Zenodo citation were corrected. The source test
suite, spelling, and URL checks also passed on the prior candidate.

The old-release check used a temporary private library containing landgraph
0.0.3; no user package library was modified. The local R-release library also
contains landgraph 0.0.3. R-hub Ubuntu R-release and macOS R-devel checks
passed on earlier candidates differing only in benchmark script packaging.
The R-hub Windows R-devel job for the prior archive flagged CRLF line endings
in its temporary Git checkout, although that archive itself contains
LF shell scripts and passed the local checks. The release branch now has a
GitHub Actions workflow. Its six jobs passed on Windows R-release and R-devel,
macOS R-release, and Ubuntu R-oldrelease, R-release, and R-devel for commit
`396790f`. Two sanitizer attempts stopped during dependency setup before
checking terradish. Local checks are not CRAN acceptance.

## Note explanation

The only NOTE is expected for a first submission: "New submission".

## Reverse dependencies

This is a first submission, so there are no CRAN reverse dependencies.
