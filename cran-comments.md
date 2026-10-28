## Resubmission

This is a resubmission of the new package terradish, first submitted as
version 1.0.0 and now version 1.1.0. In response to the CRAN review, I have:

- Added references describing the methods to the Description field, in the
  requested form: Peterman and Pope (2021) <doi:10.1111/mec.15716>, Clarke,
  Rothery and Raybould (2002) <doi:10.1198/108571102320>, and Hanks and
  Hooten (2013) <doi:10.1080/01621459.2012.724647>.
- Added `\value` sections to every exported method, including the nine
  topics listed in the review (`legacy_radish_*.Rd`,
  `legacy_terra_radish_scale_optim.Rd`, and `print.terradish_plot_list.Rd`).
  Each describes the class and structure of the result and what it means;
  the print method documents that it returns its input invisibly and is
  called for its side effect.
- Removed the `setwd()` calls that parallel workers used to load package
  sources. Workers now source files by absolute path, so no function changes
  the working directory. No package function changes `options()` or `par()`.
  The two vignette chunks that set `par()` restore it with `on.exit()`
  immediately after the change.

Since the previous submission I also added two example datasets
(`ohio_ponds`, `individual_salamanders`), four vignettes, and numerical
safeguards for nuisance-parameter fitting, with regression tests.

## R CMD check results

0 errors | 0 warnings | 1 note

- This is a new submission.

The installed size is about 8.4 Mb (data 1.3, doc 2.1, extdata 1.5, and
libs 1.6 Mb). The data and extdata hold the example datasets and the
rasters and stored results the vignettes need to run without downloads;
doc holds the twelve vignettes, whose precomputed figures keep check time
low. The help-page logo was resized to reduce the size.

## Test environments

- Windows 11 x86_64, R 4.6.1: `R CMD check --as-cran --no-manual` on the
  built archive: 0 errors, 0 warnings, and the "New submission" NOTE. A
  second, local-only NOTE reported that README.md and NEWS.md could not be
  checked because pandoc was not on the check's PATH; both convert cleanly
  with pandoc.

## Reverse dependencies

This is a new submission, so there are no reverse dependencies.
