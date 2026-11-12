# AGENTS.md

## Purpose
This repository contains an R package. The goal of any automated coding work is to make targeted, reviewable improvements while preserving package stability, documentation integrity, and reproducibility.

## General rules
- Prefer small, minimal diffs.
- Do not make broad refactors unless explicitly requested.
- Preserve existing public function interfaces unless the task explicitly requires changing them.
- Do not guess about intended statistical or scientific behavior. Flag uncertainty clearly.
- Do not silently change numerical behavior, defaults, or return values.
- When touching package code, also inspect tests, documentation, and examples that may need updating.


## Package structure
Typical files and directories that may be relevant:
- `DESCRIPTION`
- `NAMESPACE`
- `R/`
- `man/`
- `tests/testthat/`
- `vignettes/`
- `src/` and `inst/` if present
- `.github/workflows/`

## Coding expectations
- Follow the existing code style in the repository.
- Keep functions focused and avoid unnecessary abstraction.
- Add moderate annotation to code and functions for future understanding
- Prefer readable base R unless the repository already uses tidyverse-style conventions in that part of the codebase.
- Avoid introducing new dependencies unless clearly justified.
- If a new dependency is necessary, explain why and update:
  - `DESCRIPTION`
  - namespace usage
  - any affected documentation/tests

## Documentation rules
- This package uses roxygen2 documentation unless the repository clearly indicates otherwise.
- If function behavior, arguments, defaults, return values, or examples change, update the roxygen comments in `R/`.
- Do not edit generated `.Rd` files directly unless explicitly instructed.
- Keep examples lightweight and CRAN-safe.
- Do not use em-dashes within sentences. Reserve use for headings, tables, bullet lists.
- Use American English.
- Do not create or fabricate references in support of methods.

## Vignette rules
- Be as descriptive in text as possible to make the document educational.
- Provide annotation of code.
- Do not use em-dashes within sentences. Reserve use for headings, tables, bullet lists.
- Use American English.
- Do not create or fabricate references in support of methods.

## Testing rules
- For any nontrivial code change, add or update tests in `tests/testthat/`.
- Prefer focused tests that verify behavior, edge cases, and regressions.
- Do not remove failing tests merely to make checks pass unless explicitly instructed and justified.
- If behavior changes intentionally, update tests to reflect the new contract and note that clearly.

## Versioning rules
The package is released on CRAN, so the version number describes releases, not commits.

- Released versions have three parts, `MAJOR.MINOR.PATCH`:
  - `PATCH` (1.1.1 to 1.1.2): bug fixes, documentation, and performance changes that leave documented behavior and numerical results unchanged.
  - `MINOR` (1.1.1 to 1.2.0): new exported functions, arguments, models, or vignettes that are backward compatible.
  - `MAJOR` (1.1.1 to 2.0.0): changes that break existing code or change defaults, return values, or numerical results.
- Between releases, `DESCRIPTION` carries the last released version plus a fourth component, `.9000` (for example `1.1.1.9000`). Do not bump the version for individual commits.
- Choose the next release number when the release is prepared, from the changes accumulated in `NEWS.md`. Do not choose it in advance.
- The website uses this convention: a four-part version builds the development site under `dev/`, and a three-part version builds the released site.

## Updating rules
- Record every user-visible change as a brief bullet in `NEWS.md` under the top heading `terradish (development version)`. Create that heading if it is missing. Do not add a new version heading for a commit.
- Provide a succinct, minimal review of what the update contains.
- Do not commit changes until explicitly asked. Minor related edits can accumulate into one commit.

## Release rules
Apply these only when explicitly asked to prepare a release.

- Set the three-part release version and the date in `DESCRIPTION`, and rename the `NEWS.md` development heading to `terradish X.Y.Z`.
- Keep the version and date identical in `DESCRIPTION`, `CITATION.cff`, `.zenodo.json`, and `cran-comments.md`.
- Tag the released commit `vX.Y.Z` once CRAN accepts it.
- Immediately after the release, set `DESCRIPTION` to `X.Y.Z.9000` and add a new `terradish (development version)` heading to `NEWS.md`.

## Website rules
- The website is built with pkgdown from `_pkgdown.yml` and deployed by `.github/workflows/pkgdown.yaml`. Do not commit the local `docs/` preview.
- Every new exported help topic must be added to a `reference` section of `_pkgdown.yml`, and every new vignette to an `articles` section and the `Guides` navbar menu.
- Vignettes carry a two-digit reading-order number in two places, because the listings sort differently: R and RStudio sort by title, and CRAN sorts by file name. The file is `vNN-name.Rmd` and the title is `NN. Title`. The `v` is required: R rejects vignette file names that do not start with a letter. Keep the YAML `title` and `\VignetteIndexEntry{}` identical, and keep the copy in `inst/vignette-source/` in agreement with `vignettes/`.
- Avoid renumbering existing vignettes. A file name is also the `vignette("vNN-name")` name, the website URL, and the prefix of its figures in `vignettes/precomputed/`, so a renumbering must update every cross-reference in `R/`, `vignettes/`, `inst/vignette-source/`, `README.md`, `_pkgdown.yml`, and `dev/refresh-precomputed-vignettes.R`. Add new vignettes at the end of their group when possible.
- Run `pkgdown::check_pkgdown()` after changing exports, help topics, or vignettes.

## Development workflow
- Work on one feature or fix at a time, ideally on a dedicated branch.
- Make code/doc/test changes first, then inspect `git status` and `git diff` before committing.
- Run `devtools::document()`, `devtools::test()`, and `devtools::check(args = "--no-manual")` before a final commit when feasible.
- If Pandoc is unavailable in the current shell, `devtools::check(args = "--no-manual", build_args = "--no-build-vignettes", document = FALSE)` can be used as a pre-check, but run one final vignette-enabled check in RStudio when vignette files changed.
- `R CMD build` copies the whole folder, `.git` included, before applying `.Rbuildignore`. If it fails with a path-too-long error on Windows, a tool has written a deeply nested loose ref under `.git/refs/`. Run `git pack-refs --all`, which moves every ref into `.git/packed-refs` without deleting any, then remove the empty directories it leaves behind (`find .git/refs -depth -type d -empty -delete` in Git Bash; git recreates what it needs), and build again.
- After validation, make one concise commit for the batch of related changes, and merge to `master` only after that branch is stable. The version changes only at a release (see Versioning rules).

## Validation expectations
After making changes, recommend the relevant validation commands. Common commands include:

```r
devtools::document()
devtools::test()
devtools::check()
```
