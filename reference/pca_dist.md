# Pairwise genetic distances from PCA scores

Runs a PCA on allele-frequency data from an
[`adegenet::genind`](https://rdrr.io/pkg/adegenet/man/new.genind.html)
object and returns Euclidean distances computed on the retained
principal component scores.

## Usage

``` r
pca_dist(gi, n_axes = 64, scale = TRUE)
```

## Arguments

- gi:

  A `genind` object from adegenet.

- n_axes:

  Number of principal component axes to retain. Capped at the total
  number of available axes.

- scale:

  Logical. If `TRUE` (default), the resulting distance matrix is
  rescaled to `[0, 1]` via
  [`scale_to_0_1`](https://wpeterman.github.io/terradish/reference/scale_to_0_1.md).

## Value

A symmetric numeric distance matrix with one row/column per individual
in `gi`.

## Details

Missing genotypes in `gi` are imputed by the allele-mean method
(equivalent to `adegenet::tab(gi, NA.method = "mean")`), then a standard
PCA via [`prcomp`](https://rdrr.io/r/stats/prcomp.html) is applied to
the allele frequency table. Euclidean distances are then computed on the
first `n_axes` scores.

This is a convenience wrapper designed for users who prefer to work from
a `genind` object. For biallelic SNP data the more direct
[`dist_from_biallelic`](https://rdrr.io/pkg/landgraph/man/dist_from_biallelic.html)
is generally faster and does not require the adegenet package.

## See also

[`dist_from_biallelic`](https://rdrr.io/pkg/landgraph/man/dist_from_biallelic.html),
[`cov_from_genetic_data`](https://rdrr.io/pkg/landgraph/man/cov_from_genetic_data.html)

Other genetic responses:
[`check_distance_response()`](https://wpeterman.github.io/terradish/reference/check_distance_response.md),
[`lower()`](https://wpeterman.github.io/terradish/reference/lower.md)

## Examples

``` r
if (requireNamespace("adegenet", quietly = TRUE)) {
  data(nancycats, package = "adegenet")
  d <- pca_dist(nancycats, n_axes = 4)
  dim(d)
}
#> [1] 237 237
```
