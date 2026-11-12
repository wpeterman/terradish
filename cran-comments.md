## Patch release

terradish 1.1.1 is a small patch to the 1.1.0 candidate, prepared shortly
after that version. Sparse projection made graph construction impractical
on the largest benchmark rasters (above about 250,000 cells). Principal-
submatrix indexing forms the same reduced Laplacian without that bottleneck.
Numerical results and public interfaces are unchanged.

## R CMD check results

Windows 11 x86_64, R 4.6.0: `R CMD check --as-cran` on a freshly compiled
source archive, including examples, tests, vignette rebuilding, and the PDF
manual.

0 errors | 0 warnings | 1 note

The incoming-feasibility NOTE is "New submission". On November 8, 2026,
CRAN's public package index did not list terradish. The public package page
returned HTTP 404. This does not establish its manual-review queue status.

The constructor regression test passed eight expectations in 0.37 seconds.
The extended suite passed 1,814 expectations with no failures or warnings.
One source-loaded worker test was deliberately skipped in installed-package mode.
Three complete fits (melip MLPE, synthetic Wishart log-linear, and synthetic
Wishart with a Gaussian scale) reproduced 1.1.0 exactly, including standard
errors and predictions. A million-cell graph built in 14.08 seconds on one
thread. Cross-platform checks and reverse-dependency checks were not run
for this candidate. No submission to win-builder or R-hub was made.
