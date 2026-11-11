# Pairwise environmental differences from site-level covariates

Constructs pairwise environmental-mismatch covariates from site-level
values. The result can be used alongside isolation-by-resistance (IBR)
terms to model isolation by environment (IBE) within
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
or
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md).

## Usage

``` r
pairwise_endpoint_covariates(
  x,
  coords = NULL,
  transform = c("absdiff", "sqdiff", "euclidean", "manhattan"),
  scale = FALSE
)
```

## Arguments

- x:

  Site-level covariates. Supported inputs:

  - A numeric vector (one covariate per site).

  - A numeric matrix or data frame (sites in rows, covariates in
    columns).

  - A
    [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html):
    values are extracted at `coords`.

- coords:

  Required when `x` is a raster. Focal-point coordinates in the same
  projection as `x`; accepts the same inputs as
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- transform:

  How to collapse site-level values into a pairwise
  environmental-mismatch measure.

  `"absdiff"`

  : Absolute difference between endpoints: \\\|x_i - x_j\|\\. One output
    column per covariate. Symmetric; suitable when the hypothesis
    concerns an association between environmental dissimilarity and
    differentiation, regardless of direction.

  `"sqdiff"`

  : Squared difference: \\(x_i - x_j)^2\\. One output column per
    covariate. Like `"absdiff"` but penalizes large differences more
    heavily.

  `"euclidean"`

  : Euclidean distance across all covariates: \\\sqrt{\sum_k (x\_{ik} -
    x\_{jk})^2}\\. Collapses multivariate data to a single dissimilarity
    measure.

  `"manhattan"`

  : Manhattan (city-block) distance: \\\sum_k \|x\_{ik} - x\_{jk}\|\\.
    Collapses to a single measure; less sensitive to outliers than
    Euclidean.

- scale:

  Logical. If `TRUE`, standardize (z-score) site-level covariates before
  computing the pairwise transform. Recommended when covariates are on
  different scales.

## Value

A numeric matrix of class `"terradish_pairwise_covariates"` with one row
per unordered site pair and named columns reflecting the chosen
transform (e.g. `"absdiff_altitude"`, `"euclidean"`).

## Details

The function produces a \\n(n-1)/2 \times p\\ matrix of pairwise
endpoint-difference covariates, where \\n\\ is the number of focal sites
and \\p\\ is the number of output columns (one per input covariate for
`"absdiff"` and `"sqdiff"`; one for `"euclidean"` and `"manhattan"`).
Rows correspond to the lower triangle of the pairwise distance matrix,
scanned column by column.

The returned object has class `"terradish_pairwise_covariates"`. This
class carries the original site-level covariate matrix and the chosen
transform as attributes, so that
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
can automatically rebuild the correct pairwise covariates for each
train/test split.

If `x` already has class `"terradish_pairwise_covariates"`, it is
returned unchanged (pass-through).

## Interpreting pairwise effects

Construct one
[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
object to use the same environmental transforms and geographic distances
in both likelihood families. Their fitted conductance surfaces are
conditional on these terms, and conductance signs can change when the
mean structure changes. Report fits with and without pairwise effects. A
uniform `~ 1` conductance model with pairwise terms represents IBD plus
IBE, not IBE alone.

A positive environmental coefficient means differentiation increases
with environmental difference conditional on resistance. It does not
distinguish assortative mating, dispersal filtering, local adaptation,
or site variance that rises at environmental extremes. MLPE coefficients
are unconstrained; Wishart kernel weights are nonnegative. Compare their
effects in resistance-distance units with
[`terradish_ibe_ratio`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
which divides by the resistance coefficient and propagates joint
uncertainty. Ratios are undefined when that denominator is zero.

## See also

[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

## Examples

``` r
library(terra)

data(melip)
altitude <- terra::unwrap(melip.altitude)
coords <- terra::unwrap(melip.coords)

# Absolute altitude difference for every pair of sites, after
# standardizing the site values. One row per pair, in dist() order.
z_abs <- pairwise_endpoint_covariates(altitude, coords,
                                      transform = "absdiff", scale = TRUE)
head(z_abs)
#>      absdiff_spatial_template_lowres
#> [1,]                       0.5217511
#> [2,]                       1.6264202
#> [3,]                       0.6761406
#> [4,]                       0.1773603
#> [5,]                       1.1508445
#> [6,]                       0.4351686

# Several site-level columns collapsed into one Euclidean difference.
site_env <- data.frame(
  altitude = scale(as.numeric(terra::extract(altitude, coords)[, 2])),
  xcoord = scale(terra::geom(coords)[, "x"])
)
z_env <- pairwise_endpoint_covariates(site_env, transform = "euclidean")
head(z_env)
#>      euclidean
#> [1,] 0.5218581
#> [2,] 2.6250996
#> [3,] 1.5786658
#> [4,] 1.4060711
#> [5,] 1.1508931
#> [6,] 0.4478144
```
