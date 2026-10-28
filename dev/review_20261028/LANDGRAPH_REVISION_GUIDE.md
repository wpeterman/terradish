# landgraph revision guide (review of 0.0.3, 2026-10-28)

Reviewed `landgraph` on branch `main`, HEAD `49d979d`, clean working tree. Branch `codex/terradish-core-companion` sits two commits behind `main` with no unique commits, so you can delete it. Line numbers refer to that HEAD. I edited no landgraph files.

## Summary

The numerical code is sound. The Gower default, the normalized-dosage covariance, and the Hudson F_ST match their formulas; I checked F_ST against a hand calculation. The help text and the vignette fall short. They contradict terradish on `nu`, describe a `diagonal` default that 0.0.3 removed, and never explain why terradish divides `cov_from_biallelic()` by `N[1]`. terradish re-exports these topics, so its users meet the contradictions under terradish's own name. Three small code defects affect what those users see: `dist_from_cov()` drops row names, `fst_from_biallelic()` drops names and silently misaligns a locus-length `N`, and `dist_from_biallelic()` requires `N`. Several pages also describe directed models that terradish does not provide and claim that the two packages' graph objects are interchangeable, which they are not. Ship the fixes as landgraph 0.0.4, but do not hold terradish 1.0.0 for that release.

## Critical

### ?cov_from_genetic_data: `nu` (R/genetic_covariance.R:306-320)

Current: "`nu` is usually the retained SNP count ... use the number of loci". Replace the paragraph with:

```
Wishart models in terradish, \code{terradish::wishart_covariance()} and
\code{terradish::generalized_wishart()}, need \code{nu}, an
effective-information value that you supply. This function does not compute
it, and the number of SNPs, loci, or allele features does not supply it:
linked loci share population history, so they carry far less independent
information than their count suggests. Hold \code{nu} fixed across compared
models and report sensitivity with \code{terradish::terradish_rescale_nu()}.
```

Rationale: the current text contradicts every terradish guide and documents an argument this function does not have.

### Vignette: stale defaults and `nu` (vignettes/getting-started.Rmd)

- 261-266, current "defaults to `diagonal = "within"`" (the knitted output below it prints `"gower"`). Replace with: "Because `cov_from_genetic_data()` defaults to `diagonal = "gower"`, the diagonal holds each population's variance on the same scale as the covariances, rows sum to zero, and the matrix is positive semidefinite, as a Wishart likelihood requires. Population-graph workflows instead use `diagonal = "within"`, whose diagonal is on another scale and does not suit Wishart models."
- 288-290, current "`"auto"` (the default)". Replace with: "`"gower"` (the default) keeps the double-centered diagonal; `"within"` substitutes within-population variance; `"auto"` restores the default used before 0.0.3."
- 274 and 282: state that `within_variance` fills the diagonal only with `diagonal = "within"`.
- 326-332, current "For SNPs that is roughly the retained SNP count". Replace with the paragraph above, in Markdown.

Rationale: the prose contradicts the code and the output printed beneath it.

### Vignette: unclosed code fence (getting-started.Rmd:485, 506-507)

Line 485 ends with "`# symmetric resistance```", so the chunk never closes. In my knitted copy, the summary table, See also, and References all render as code. Put the fence on its own line, and add a blank line before `# References`.

## Major

### ?dist_from_cov (R/genetic_covariance.R:656-679, 698-702)

- Code: the result has column names but no row names, which is why terradish's getting-started output prints `[1,]` beside `P01`. It also inherits every attribute of `Cov`. See the code section.
- 656-658, current "inverts the Gower double-centering". Replace with: "Converts a covariance to the squared Euclidean distances it implies. For \code{cov_from_genetic_data()} output with \code{diagonal = "gower"}, the result equals its \code{centroid_distance2} attribute; with \code{"within"}, it does not."
- 673-676: "Squared distances from a coherent covariance suit \code{terradish::generalized_wishart()}, which gives the same likelihood as \code{terradish::wishart_covariance()} and needs the same \code{nu}. \code{terradish::mlpe()} accepts any pairwise distance and needs no \code{nu}."
- 678-679, current "non-negative off-diagonal entries". Replace with: "A symmetric matrix with the dimnames of \code{Cov} and a zero diagonal. Off-diagonal entries are nonnegative when \code{Cov} is positive semidefinite."

### ?cov_from_biallelic (R/genetic_covariance.R:73-74, 113-118, 130)

Current (73-74): "directly suitable as the response \code{S} in \code{wishart_covariance}". Delete that sentence and add to Details:

```
The function standardizes counts, not frequencies, so entries grow with the
haploid sample size: \eqn{Z_{i\ell} = \sqrt{N_{i\ell}}(\hat p_{i\ell} -
p_\ell)/\sqrt{p_\ell(1-p_\ell)}}. For individual 0/1/2 dosage this is the
genomic relationship matrix scale. For population counts with equal sample
sizes, \code{cov_from_biallelic(Y, N) / N[1]} gives the standardized
allele-frequency covariance that terradish fits. Do not use that shortcut
with unequal sizes; rows then carry unequal sampling variance that a common
nugget does not absorb. \code{Y} and \code{N} cannot contain \code{NA}. For
population counts, record missing genotypes by lowering the matching cells
of \code{N}. Every cell of \code{N} must be positive, so for individual
dosage you must filter loci or individuals first. Counting either allele
gives the same covariance.
```

At line 115, after "report the filtering rule", add "and check that conclusions hold across thresholds". At line 130, describe the result as an `nrow(Y)` by `nrow(Y)` matrix on the count-standardized scale.

Rationale: terradish divides by `N[1]` and warns against that shortcut with unequal sizes, but this page never says why.

### ?cov_from_genetic_data: missing data (R/genetic_covariance.R:300-304 vs 333-334)

Current (333-334): "Missing values are not currently supported", yet lines 302-304 describe how missing allele calls are imputed. Replace 333-334 with:

```
A numeric feature matrix cannot contain \code{NA}. For allele calls, each
missing copy becomes the most frequent allele at its locus; the
\code{imputed_allele_calls} attribute records the count. Modal imputation
shrinks distances, so report it. \code{terradish::pca_dist()} imputes allele
means instead. Code missing calls as \code{NA}: the function reads codes such
as 0 or -9 as real alleles.
```

Rationale: in my test, a call coded 0 became the feature `l1:0` with no warning.

### ?cov_from_genetic_data: arguments and \value (R/genetic_covariance.R:280-297, 357-383)

- `center`: "The covariance depends only on differences among units, so centering changes only the \code{centroids}, \code{unit_features}, and \code{feature_center} attributes." I verified that the covariance is unchanged to within 7e-16.
- `scale`: "Divides each feature by its standard deviation, so rare alleles get the same variance as common ones and therefore more relative weight. Constant features are dropped with a warning."
- `diagonal`: "\code{"gower"} (default) keeps the double-centered diagonal, so every entry shares one scale; use it for Wishart models. \code{"within"} substitutes summed within-population variances, as in population graphs; terradish warns when it receives one. \code{"auto"} restores the default used before 0.0.3."
- `normalize`: append "The divisor counts expanded allele features, not loci, and is not \code{nu}."
- \value 367-369 is wrong: `within_variance` is attached whenever `groups` is given, with 0 for single-member groups. Lines 370-371: `unit_size` is numeric for groups. Also document `unit_features`, `unit_distance2`, `input`, `normalize`, and the `imputed_*` attributes, and add "Use \code{S[, ]} to view the matrix without attributes."

### ?dist_from_biallelic (R/genetic_covariance.R:710-718, 732-735)

- Code: `N` has no default, so `dist_from_biallelic(Y)` fails, even though the help points to the flexible `N` forms. Change the signature to `(Y, N = NULL, ...)` and pass `...` to `cov_from_biallelic()`.
- 716-718: "Squared distances on the count scale of \code{cov_from_biallelic()}. With equal sample sizes, terradish's response \code{dist_from_cov(cov_from_biallelic(Y, N) / N[1])} equals this matrix divided by \code{N[1]}."

Rationale: users following terradish otherwise get distances `N[1]` times larger than the vignette's `D`.

### ?fst_from_biallelic (R/genetic_covariance.R:9-11, 32-33, 46-66)

- Code: the result has no dimnames. A scalar `N` works by accident. A locus-length vector recycles down the columns with no warning; in my test, the result was 44% off. Route `N` through `.biallelic_sample_size_matrix()` and copy the row names.
- `@param N`: "Haploid sample sizes in any form accepted by \code{cov_from_biallelic()}; every entry must exceed 1."
- `@return`: "A symmetric matrix named by \code{rownames(Y)} holding pairwise Hudson F_ST, with a zero diagonal."
- Lines 22-25 already match terradish (use `mlpe`; not a Wishart response). Keep them.

### ?deme_graph, ?edge_gradient, vignette, README: terradish claims

- R/deme_graph.R:5-7, vignette 90-92, README 23-25, current "interchangeable with a \code{terradish_graph}". Replace with: "The class lets \code{edge_gradient()} and \code{edge_flow()} accept the object. It lacks the covariates, Laplacian, and Cholesky factor that \code{terradish::terradish()} needs; build fitting graphs with \code{terradish::conductance_surface()}." I did not run `terradish()` on a `deme_graph`. The missing fields come from reading terradish R/radish_graph.R:284-297.
- R/edge_covariates.R:4 and 15 cite `terradish_directed`, which does not exist. Lines 6-7 describe a coefficient in a model that terradish does not fit. Replace 3-7 with: "Builds the antisymmetric covariate \eqn{d_{ab} = x_a - x_b} for directional gene-flow models; it is positive from high to low \code{x}." Fix vignette 394-395 and README 54-57 the same way.

## Minor

- American English: R/deme_graph.R:13, 15, 17 "neighbour"; vignette 111 and 453 "grey", 115 "neighbourhood", 136 "neighbour". Keep the argument `neighbours`, and prune inst/WORDLIST to match.
- DESCRIPTION:11-18: cite the methods the Rd files already reference: Yang et al. (2010) <doi:10.1038/ng.608>, Gower (1966) <doi:10.1093/biomet/53.3-4.325>, Dyer and Nason (2004) <doi:10.1111/j.1365-294X.2004.02177.x>, Bhatia et al. (2013) <doi:10.1101/gr.154831.113>. All four were verified in Consensus on 2026-10-28 (authors, year, title, and DOI match). Drop "No compiled code." and update `Date`.
- inst/CITATION:7 hard-codes "0.0.1". Use `sprintf("R package version %s", meta$Version)`.
- R/deme_graph.R:19-22 omits the `coords` element. Lines 81-82 do not validate `k`: `k = 0` returns three spurious edges. Line 14 says "integer grid", but the code needs only uniform spacing.
- Write bare references to `wishart_covariance`, `generalized_wishart`, `mlpe`, and `simulate_covariance_response` (lines 74, 140, 320, 385-386, 675-676, 717-718) as `\code{terradish::name()}`, without `\link`.
- R/edge_covariates.R:12-14: state that rows m+1 to 2m reverse `edge_pairs`.
- Vignette 255: use `round(Sg[, ], 3)`, because the current call prints 15 attributes. Line 198: replace "exactly the response S" with the `S / 40` conversion. Line 211: update after the `fst` fix.
- Optional, line 60: `(outer(r, r, "+") - 2 * tcrossprod(Fr)) / L` with `r = rowSums(Fr^2)` gives the same `f2` without the n^2 by L array (verified).
- README 14-17: install from CRAN, and define `Y` and `wind_field` in the example.

## Code and test changes

```r
dist_from_cov <- function(Cov)
{
  Cov <- as.matrix(Cov)
  if (!is.numeric(Cov) || nrow(Cov) != ncol(Cov))
    stop("`Cov` must be a square numeric matrix", call. = FALSE)
  d <- diag(Cov)
  D <- outer(d, d, "+") - 2 * Cov[, , drop = FALSE]  # drops other attributes
  dimnames(D) <- dimnames(Cov)
  # terradish warns when attr(S, "diagonal") is "within"; keep that label.
  attr(D, "diagonal") <- attr(Cov, "diagonal", exact = TRUE)
  attr(D, "centered") <- attr(Cov, "centered", exact = TRUE)
  D
}
# fst_from_biallelic(): N <- .biallelic_sample_size_matrix(N, Y, ploidy = 2)
#   and dimnames(fst) <- list(rownames(Y), rownames(Y))
```

Keeping only `diagonal` and `centered` is a design choice for you to confirm. Add tests/testthat/test-distance-fst.R. I ran these tests against patched copies: all pass, and the final `fst` expectation fails on the current code.

```r
test_that("dist_from_cov keeps names and inverts Gower", {
  x <- matrix(c(0, 1, 1, 1, 2, 0, 2, 1, 0, 2, 1, 2), ncol = 2, byrow = TRUE)
  C <- cov_from_genetic_data(x, groups = rep(c("p1", "p2", "p3"), each = 2))
  D <- dist_from_cov(C)
  expect_identical(dimnames(D), dimnames(C))
  expect_equal(D[, ], attr(C, "centroid_distance2"))
  Y <- matrix(c(2, 1, 0, 1, 1, 1, 0, 1, 2), 3, byrow = TRUE)
  expect_equal(dist_from_biallelic(Y), dist_from_cov(cov_from_biallelic(Y)))
})

test_that("fst_from_biallelic matches Hudson and aligns N", {
  set.seed(1)
  Y <- matrix(rbinom(15, 16, 0.4), 3, dimnames = list(c("p1", "p2", "p3"), NULL))
  N <- matrix(20, 3, 5); p <- Y / N
  num <- (p[1, ] - p[2, ])^2 - p[1, ] * (1 - p[1, ]) / 19 - p[2, ] * (1 - p[2, ]) / 19
  den <- p[1, ] * (1 - p[2, ]) + p[2, ] * (1 - p[1, ])
  fst <- fst_from_biallelic(Y, N)
  expect_equal(fst[1, 2], mean(num) / mean(den))
  expect_identical(rownames(fst), rownames(Y))
  expect_equal(fst_from_biallelic(N - Y, N), fst)
  Nl <- c(20, 20, 18, 20, 16)
  expect_equal(fst_from_biallelic(Y, Nl),
               fst_from_biallelic(Y, matrix(Nl, 3, 5, byrow = TRUE)))
})
```

Also add a test that `center = TRUE` and `center = FALSE` give equal covariances, and that `deme_graph(..., "knn", k = 0)` errors.

## Release and coordination

1. Branch from `main`. Bump to 0.0.4 and update `Date`. NEWS: the documentation now matches terradish (`nu`, `N[1]` scaling, the `diagonal` default, missing data). `dist_from_cov()` keeps dimnames. `fst_from_biallelic()` keeps names and accepts flexible `N`. `dist_from_biallelic()` gains `N = NULL`. `deme_graph()` validates `k`. Values change only for the formerly misaligned locus-length `N`.
2. Run `devtools::document()`, `devtools::test()`, `devtools::spell_check()`, `R CMD check --as-cran` on the built tarball, and `devtools::check_win_devel()`.
3. Install the patched landgraph and run terradish's tests and check, because terradish calls `dist_from_cov()` internally and in several tests.
4. Submit terradish 1.0.0 first. Its `landgraph (>= 0.0.3)` stays correct: doc fixes need no higher minimum, since installation pulls the current CRAN landgraph. Until 0.0.4 lands, the re-exported help still shows the old `nu` text.
5. Submit landgraph 0.0.4 after CRAN accepts terradish. CRAN asks for updates no more than every one to two months, and 0.0.3 is dated 2026-09-24, so aim for late November 2026 or later. List terradish as a checked reverse dependency in cran-comments.md.
6. Raise terradish to `landgraph (>= 0.0.4)` only when its code, tests, or precomputed output rely on the new names or the `N` default, then refresh the precomputed getting-started output.

## Checklist

- [ ] Replace the `nu` text in ?cov_from_genetic_data and the vignette
- [ ] Fix the vignette `diagonal` default and the unclosed fence
- [ ] Document count scale, `N[1]`, and missing data in ?cov_from_biallelic
- [ ] Reconcile missing data and `NA` coding; rewrite four arguments and `\value` in ?cov_from_genetic_data
- [ ] Fix `dist_from_cov()` names and help
- [ ] Add `N = NULL` to `dist_from_biallelic()` and state its scale
- [ ] Fix `fst_from_biallelic()` names and `N`
- [ ] Remove the `terradish_directed` and interchangeability claims
- [ ] Clear the minor items
- [ ] Add tests; run document, test, spell check, and `--as-cran`
- [ ] Run terradish against the patched landgraph
- [ ] Bump to 0.0.4; update NEWS and cran-comments
- [ ] Submit terradish 1.0.0, then landgraph 0.0.4
