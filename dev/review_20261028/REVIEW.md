# terradish review: `codex/amphibian-vignette`, October 28, 2026

Scope: the uncommitted working tree in
`C:/Users/peterman.73/.codex/worktrees/amphibian-vignette/terradish`
(Codex's completed integration, version 1.0.0). Standards: the
peterman-research-suite skills `r-package-docs`, `r-doc-style`,
`r-function-documentation`, `r-vignette-writing`, and `r-cran-compliance`.
Four parallel reviewers read all 11 vignettes and the two dataset help
topics; every finding acted on below was re-verified against the source,
the printed output, or the figures. Nothing is committed.

A pre-review backup of the tree (excluding `.git` and `dev/`) is in this
session's scratchpad as `worktree_backup_before_review.tgz`.

## 1. CRAN requests

| Request | Status | Evidence |
|---|---|---|
| Method references in DESCRIPTION | Done, extended | Clarke, Rothery and Raybould (2002), Hanks and Hooten (2013) were already present. Added Peterman and Pope (2021) <doi:10.1111/mec.15716>, which describes the core method (one parameterized conductance surface linked to genetic data). All three verified in Consensus. Consensus lists the online year (2020); the issue year 2021 is used. |
| `\value` on the nine flagged Rd files | Done | All nine now describe class, structure, and meaning (Codex). I also expanded `terradish()`, `terradish_grid()`, and `terradish_distance()`, whose `\value` named a class without saying what it meant ("An object of class terradish_grid"). Every exported topic now has `\value`. |
| `options()`, `par()`, `setwd()` | Done | The two `setwd()` calls in `R/radish_grid.R` (HEAD) are gone; workers source files by absolute path. No package function calls `par()` or `options()`. The two vignettes that set `par()` restore it with `on.exit()` inside `local()`. Tests that change the working directory restore it immediately. |

`cran-comments.md` was an internal Codex log (it named the branch and
local receipts). It is rewritten as a resubmission note listing the three
fixes; the original is preserved at
`dev/vignette_consistency_20261028/cran-comments-development-record.md`.
Fill in the test-environment lines from the final check.

### Other CRAN-facing problems found and fixed

- **URL check failure:** `https://doi.org/10.1111/mec.12747` returns 403
  (Wiley blocks checkers) and would be reported from the built vignettes.
  Now plain `doi:` text.
- **Encoding defect:** `large-landscapes.Rmd` contained a double-encoded en
  dash that rendered as "15â€“23%". Fixed.
- **Shipped internal notes:** `inst/extdata/individual-salamanders/provenance/`
  contained a compute-budget and "October 23 logoff" paragraph and 40
  absolute paths under `C:\Users\peterman.73\.codex\...`. Removed and made
  relative. The rest of `SAMPLING_PROTOCOL.md` still reads as an internal
  work order; consider rewriting it as a plain sampling description.
- **Duplicate files:** three unreferenced PNGs in `vignettes/precomputed`
  duplicated the extdata copies; removed.

### Before you commit or build in RStudio

- `dev/` in this worktree is **33 GB**, almost all untracked Codex
  artifacts (ten nested `terradish.Rcheck` folders, receipts, libraries).
  Do not `git add dev/` wholesale; stage only the files you want tracked.
- `R CMD build` (and so `devtools::build()`/`check()`) copies the whole
  folder before applying `.Rbuildignore`. The nested `dev/**/terradish.Rcheck`
  paths exceed the Windows 260-character limit, so building directly from
  this worktree fails with "copying to build directory failed". Delete or
  move those check folders, or build from a copy without `dev/` (as done
  for the validation below).

### Still open for CRAN

- **Size.** The checked tarball doubled from 2.9 MB (HEAD) to 6.4 MB, and the
  installed size is 9.4 MB (data 1.3, doc 1.9, extdata 1.5, help 1.6,
  libs 1.6 MB). CRAN's Linux checks will likely report an installed-size
  NOTE. Either trim (the 0.9 MB `predictions.rds`, precomputed figure
  resolution) or explain it in `cran-comments.md`.
- **Version.** Kept at 1.0.0 because the earlier submission was not
  accepted; AGENTS.md asks for a bump at commit time.

## 2. Package code change

`plot(type = "marginal")` and `plot(type = "marginal_response")` set the y
title through `scale_y_continuous(name = ...)`, which silently overrides a
user's `labs(y = ...)`. Every vignette that relabeled the y axis showed the
wrong label in its figure. The title now comes from `labs()`. Three tests
were updated and a regression test added (`test-s3-methods.R`). Default
labels are unchanged.

## 3. Function documentation

- `?terradish`: `@return` now lists every element (`convergence`, `fit`
  with `phi`, `fitted`, `subproblem`, `loglik`, `mle`, `comparison`, and
  the rest). It previously pointed to a nonexistent `fit$subproblem`, and
  the vignettes used undocumented internals.
- `?mlpe`: said `summary()` reports the shared-site correlation on its
  natural scale. It prints logit `rho`. Corrected, with the conversion
  `plogis(rho) / 2` and readings for `alpha`, `beta`, and `tau`.
- `?wishart_covariance`: removed a changelog paragraph ("The former
  full-rank likelihood ... has been removed").
- `?plot.terradish`: notes that `marginal_response` with a single covariate
  is the conductance curve inverted and rescaled.
- `?terradish` details: removed a `--` dash pair (em-dash rule).

**Outside this package:** landgraph's `?cov_from_genetic_data` says `nu`
"is usually the retained SNP count," which contradicts every terradish
guide. Fix it in landgraph.

## 4. Vignettes

The guides are accurate and well organized, but they share three
weaknesses: stacked defensive caveats ("this is not X, not Y"),
development-history voice, and printed output that is never read in
numbers. Reviewers estimated 20 to 40% of the prose could go. This pass
corrected every error, added a "How to read" reading with numbers for
each key output, cut the clearly redundant passages, and completed the
cross-links. A deeper concision pass remains (see section 6).

The seven precomputed guides were re-rendered from their edited executable
sources against the package build in this tree (installed to a scratch
library, not your R library). Every fitted coefficient, log-likelihood,
and cross-validation score is identical to Codex's receipts (maximum
difference 0), and all 23 fits converged. Prose numbers are computed
inline in the sources, so they cannot drift from the output. Total length
is roughly unchanged: cuts were offset by the new readings.

| Guide | Main corrections |
|---|---|
| getting-started | Quick reference failed (`ggplot2` not loaded, graph built with different arguments); now runnable and complete. Added pairwise-distance (MLPE) pointer, plain `nu` explanation, SD in forest units, ratio interval (`exp(confint())`), AIC decision rule, `converged` column explanation. Removed duplicated figure block with one caption repeated under four images, a CV paragraph in a guide with no CV, the "reading times" sentence, `melip` changelog line, and an R-devel link. |
| model-comparison | **The LRT reading said the p-value was "small"; it is 0.78.** Rewritten. Explained `mean_difference`, `se_difference`, and `within_one_se`, including why `forest` is flagged `FALSE` while `full` is `TRUE`. AIC/AICc/BIC decision rule with numbers. CV plot lines now connect candidates in order of complexity. Quick reference now includes `forest_slope`. Typos fixed. |
| gaussian-scale-optimization | Added a reading of `gaussian_scale_summary()` (sigma in meters and cells, radial 50/90% distances). Kernel caption now says the band is a Wald interval, not the profile interval. Profile, AIC, and CV readings with decision rules. Removed unverifiable development history ("347 and 518 m"). |
| spline-conductance | Genetic-response reading called a covariance plot "genetic distance"; fixed. Read the fitted curve, focal support, and `cv$summary`. Corrected the basis guidance (`s(x, df, basis = "bs")`). Changelog voice removed. |
| categorical-conductance | Estimate now compared with the known true ratio (0.37; interval 0.06 to 0.71 contains it). CV verdict stated (24 ponds cannot distinguish the class effect from uniform). Off-topic data-generation opener replaced. |
| ibe-ibr-workflow | IBR and IBE defined. **Added the only full MLPE summary reading in the package** (`alpha`, `beta`, environmental slope, `tau` as residual SD, `rho` as shared-pond correlation 0.47) and the resistance-distance equivalent table (`terradish_ibe_ratio()`). AIC tables, `mlpe_response_change()`, and CV read in numbers. |
| simulation-design | Fit plot reading corrected (least-squares line, variance and covariance clusters). Recovery stated against the generating 0.6. Startup messages suppressed. |
| wishart-covariance | Added the missing decision rule: covariance (Wishart) versus pairwise distance (MLPE). Heatmap no longer washed out by the diagonal. |
| large-landscapes | Encoding defect fixed. Added a demonstration of `terradish_assess_settings()`, `terradish_solver_benchmark()`, and `crop_to_focal_buffer()` (eval = FALSE) with a reading. Stranded paragraph removed. |
| individual-salamanders | Parameter table showed code that errors (`phi["tau", "Estimate"]`). Quick reference used undefined objects; now self-contained. Variable shadowing `conductance()` renamed. IBE section that fits nothing condensed. Peak formula explained. Development-history text cut. |
| data-generation | Development-log passages cut, undefined jargon (non-Wright-Fisher, recapitation, JC69) explained, the counts actually read. |

All 11 guides now link back to getting-started and forward to the guides
a reader needs next.

## 5. Critical coverage gap: the pairwise-distance MLPE workflow

The CRAN-submitted `getting-started` taught the standard workflow most
users arrive with: a pairwise genetic distance (FST, Dps, individual
distances) fitted with `measurement_model = mlpe`. The rewrite replaced it
with a Wishart covariance built from genotypes. After this pass, MLPE is
covered by a pointer in getting-started, a decision rule in
wishart-covariance, and a full summary reading in ibe-ibr-workflow. That
reading is embedded in an IBE negative control, though. No guide yet walks
a distance-matrix user from their matrix to an interpreted MLPE fit. A
short dedicated guide (the `melip` data are ideal and fast) is the
highest-value addition remaining.

## 6. Recommended next steps (not done)

1. Add the MLPE pairwise-distance guide (section 5).
2. Second concision pass, about 25% of prose, mainly in getting-started
   (the "Meet the forest-pond example" section duplicates data-generation),
   ibe-ibr-workflow (the optional kernel block), large-landscapes (the `nu`
   section duplicates wishart-covariance), and individual-salamanders.
3. **Scientific check:** in getting-started the simulated forest effect is
   6 per unit forest, or 0.79 per SD, but the fitted coefficient is 0.58
   with a 95% interval of 0.40 to 0.76 under `nu = 40`. The truth sits just
   outside the interval. The generating exchange model differs from the
   fitted Wishart model, so exact recovery is not expected, but decide
   whether and how to say so.
4. gaussian-scale-optimization starts the fit at the true 500 m; start at
   300 m and re-render.
5. individual-salamanders: the worked held-out fold (fold 1) gains mostly
   from matching the mean distance (r = 0.02 for both models); fold 2 shows
   real pair-level gains. The three figures are static PNGs whose plotting
   script is not shipped.
6. simulation-design shows no `covariance_response_power()` output; a small
   precomputed sweep would teach its columns.
7. `?ohio_ponds` never mentions `extdata/ohio-ponds/landcover100.rds` (the
   categorical guide's only raster), and the `simulation` fields are not
   listed. `?individual_salamanders` leaves several record fields undefined.

## 6b. Second round (same day): MLPE guide, concision, size, landgraph

**New guide: `vignette("mlpe-pairwise-distance")`.** Closes the coverage
gap in section 5. It uses all 37 `melip` populations at full resolution
(22,443 cells), is precomputed like the other computational guides (about
one minute to render), and teaches: when to use MLPE; why
`check_distance_response()` rejects an FST matrix with negative entries;
sites that share a raster cell; IBD, forest, and forest + altitude fits;
the full MLPE summary (conductance ratio 1.74 per SD, 95% CI 1.47 to
2.07; `beta`; `alpha`; residual SD; shared-site correlation 0.29); how
forest's coefficient collapses (0.56 to 0.08) once altitude enters;
the likelihood-ratio test; and spatial cross-validation, where the
altitude model's 74-unit AIC advantage fails to transfer (it loses 28.7
units in one block). Links to it were added from getting-started,
wishart-covariance, model-comparison, ibe-ibr-workflow,
individual-salamanders, `?mlpe`, `?melip`, `?terradish-package`, and the
README.

**Concision pass.** Three agents edited the guides under hard
constraints; an independent check confirmed every code chunk, inline `r`
expression, and YAML block is unchanged.

| Guide | Prose words before | After | Cut |
|---|---|---|---|
| getting-started | 2,482 | 2,161 | 12.9% |
| gaussian-scale-optimization | 1,744 | 1,532 | 12.2% |
| spline-conductance | 1,358 | 1,189 | 12.4% |
| large-landscapes | 1,951 | 1,720 | 11.8% |
| individual-salamanders | 2,924 | 2,671 | 8.7% |
| categorical-conductance | 1,106 | 1,013 | 8.4% |
| model-comparison | 1,520 | 1,402 | 7.8% |
| data-generation | 1,502 | 1,402 | 6.7% |
| ibe-ibr-workflow | 1,665 | 1,564 | 6.1% |
| simulation-design | 931 | 888 | 4.6% |
| wishart-covariance | 983 | 945 | 3.9% |
| **Total** | **21,090** | **19,158** | **9.2%** |

The 20 to 30% target was not reached without deleting readings,
decision rules, or one-off scope caveats. Remaining candidates for deeper
cuts, if you accept losing them: the data-generation equilibrium, RADseq,
and second-assay caveats; getting-started's "adding a constant to S" and
"values are not probabilities of movement" sentences; the pedigree-check
paragraphs in gaussian-scale-optimization; the comparison-contract
paragraph in model-comparison.

**Package size.** Done without removing content: the help logo
(`man/figures/terradish-sticker.png`, 956 x 1120 px) was resized to
400 px wide, 1.34 MB to 76 KB (the original is in git history and the
session backup); three duplicate PNGs were removed earlier. Measured but
not applied, because each removes or alters data:

| Option | Saves | Cost |
|---|---|---|
| Drop `ohio_ponds$raw_snp` and `raw_microsat` (exact supersets of `snp` and `microsat`) | 273 KB | Raw-capture genotypes no longer shipped; update `?ohio_ponds` and its test |
| Round `predictions.rds` values to 6 significant digits | 408 KB | Maximum change 5e-6; lossy |
| Keep only the 2 of 40 prediction frames the guide reads | 862 KB | Other held-out predictions no longer shipped and not regenerable from shipped code |
| Palette-reduce figure PNGs | about 120 KB | Slight color banding |
| Lossless xz recompression of extdata `.rds` | 0 | Already xz |

The installed size will stay above CRAN's 5 MB NOTE threshold even with
all of these (libs, help, doc, and data each exceed 1 MB), so explain it
in `cran-comments.md` rather than cutting content.

**landgraph.** `dev/review_20261028/LANDGRAPH_REVISION_GUIDE.md` gives
prioritized revisions with replacement text, code patches, tests, and a
release order (terradish 1.0.0 first, landgraph 0.0.4 after acceptance).
Critical items: `?cov_from_genetic_data` tells users `nu` is the SNP count
(verified), the landgraph vignette has an unclosed code fence (verified:
odd fence count) and a stale `diagonal` default, `dist_from_cov()` drops
row names, and `fst_from_biallelic()` misaligns a per-locus `N`. The four
DESCRIPTION citations it proposes were verified in Consensus.

## 7. Validation

- `devtools::document()` (roxygen2 8.0.0): regenerated `terradish.Rd`,
  `terradish_grid.Rd`, `terradish_distance.Rd`, `mlpe.Rd`,
  `wishart_covariance.Rd`, and `legacy_radish_fit_wrapper.Rd`; all pass
  `tools::checkRd()`. `plot.terradish.Rd` was regenerated after the
  package check below (a one-paragraph help-text addition only).
- Touched test files with `NOT_CRAN=true` (gaussian-scale-conductance,
  s3-methods, session-state, wishart-covariance): 25 tests, 201
  expectations, 0 failures, errors, or skips.
- Seven precomputed guides re-rendered; results identical to the earlier
  receipts (section 4).
- `R CMD build` from a copy without `dev/`, then
  `R CMD check --as-cran --no-manual`, Windows 11, R 4.6.1:
  **0 errors, 0 warnings, 2 NOTEs**. The NOTEs are "New submission" and a
  local-only note that README.md and NEWS.md could not be checked without
  pandoc on the PATH (both convert cleanly with pandoc). Remote URL checks
  reported no problems. Installed size 9.4 Mb (INFO). Archive 6.27 MB.
  Timings: install 80 s, examples 20 s, `--run-donttest` examples 113 s,
  tests 49 s, vignette rebuild 38 s.
- Not run: R-devel, macOS, Linux, win-builder, PDF manual. Run
  `devtools::check_win_devel()` before submitting.

**Second round validation** (after the MLPE guide, concision pass, and
size work):

- All eight precomputed guides re-rendered; every fit, log-likelihood,
  and CV score is identical to the earlier receipts, and the three new MLPE
  fits converged in both optimizers.
- Mechanics: 12 guides, titles match their index entries, every figure
  exists, no broken links or non-ASCII characters.
- `devtools::document()` regenerated `individual_salamanders.Rd`,
  `melip.Rd`, `mlpe.Rd`, and `terradish-package.Rd`.
- `R CMD check --as-cran --no-manual`, Windows 11, R 4.6.1: **0 errors,
  0 warnings, 2 NOTEs** (New submission; local-only pandoc note).
  Archive 5.19 MB (was 6.27 MB). Installed size 8.4 Mb (was 9.4 Mb);
  help is no longer above 1 Mb. Timings: install 81 s, examples 22 s,
  `--run-donttest` 143 s, tests 62 s, vignette rebuild 43 s.
