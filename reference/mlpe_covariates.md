# MLPE measurement model with pairwise environmental covariates

Creates an MLPE measurement model whose mean structure includes both an
IBR term (resistance distance) and one or more fixed pairwise
environmental covariates (IBE terms).

## Usage

``` r
mlpe_covariates(
  x,
  coords = NULL,
  transform = c("absdiff", "sqdiff", "euclidean", "manhattan"),
  scale = FALSE,
  rho_start = 0.2
)
```

## Arguments

- x:

  Pairwise endpoint-difference covariates (the output of
  [`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md)),
  *or* site-level environmental data in any format accepted by that
  helper (numeric vector, matrix, data frame, or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)).

- coords:

  Optional focal-point coordinates, required when `x` is a raster.
  Accepts the same formats as
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- transform:

  Pairwise transform applied when `x` is not already a
  `"terradish_pairwise_covariates"` object. Passed through to
  [`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md).

- scale:

  Logical. Standardize site-level covariates before the pairwise
  transform? Passed to
  [`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md).

- rho_start:

  Starting value for the MLPE correlation parameter \\\rho\\ on its
  natural scale, which must lie in `(0, 0.5)`. The default `0.2`
  corresponds to a moderate shared-site correlation. Adjust if the
  optimizer fails to improve from the starting value.

## Value

A function of class `"terradish_measurement_model"` suitable for use as
the `measurement_model` argument of
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
and
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md).

## Details

The fitted mean structure is:

\$\$S\_{ij} = \alpha + \beta R\_{ij} + Z\_{ij}\gamma + e\_{ij}\$\$

where \\R\_{ij}\\ is the resistance distance (IBR) and \\Z\_{ij}\\ is
the row of pairwise endpoint-difference covariates for the pair
\\(i,j)\\ (IBE). The residuals follow the MLPE correlation structure of
Clarke et al. (2002).

The \\\gamma\\ coefficients are conditional on resistance distance, the
other pairwise covariates, and the MLPE correlation model. They describe
associations with the genetic-distance response. They do not cleanly
separate IBE from IBR, establish causality, or identify a unique
ecological mechanism when spatial predictors are correlated.

The nuisance parameters (intercept \\\alpha\\, IBR slope \\\beta\\, IBE
coefficients \\\gamma\\, log-precision `tau`, and the unconstrained
correlation parameter `rho`) are all profiled out during optimization;
only the conductance parameters in the formula are optimized by the
outer loop. The actual MLPE shared-site correlation is
\\\text{plogis}(\rho)/2\\.

When `x` is site-level data or the output of
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
the pairwise covariate matrix is stored as an attribute of the returned
function. This allows
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
to rebuild the correct pairwise covariates for each train/test split
automatically.

Unclassed matrices are interpreted as pairwise columns and must have
\\n(n-1)/2\\ rows. Use a data frame for site-level columns, or
[`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md)
to construct their differences.
[`pairwise_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md)
combines these with other pairwise inputs.

## Worked example

See
[`vignette("v08-ibe-ibr-workflow", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v08-ibe-ibr-workflow.md)
for paired environmental differences, genetic-response contrasts, and
their distinction from covariance kernels.

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

## References

Clarke RT, Rothery P, Raybould AF (2002). Confidence limits for
regression relationships between distance matrices: estimating gene flow
with distance. Journal of Agricultural, Biological, and Environmental
Statistics, 7, 361–372.
[doi:10.1198/108571102320](https://doi.org/10.1198/108571102320) .

## See also

[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

## Examples

``` r
# \donttest{
library(terra)

data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)
names(covariates) <- c("altitude", "forestcover")
# Coarsened so the example runs quickly; name layers, then scale.
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
coords <- coords[keep]
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords)

# Site-level altitude becomes a pairwise absolute difference (IBE term).
g_joint <- mlpe_covariates(terra::unwrap(melip.altitude), coords,
                           transform = "absdiff", scale = TRUE)
fit_joint <- terradish(Fst ~ altitude + forestcover, surface,
                       measurement_model = g_joint)
# The nuisance table adds absdiff_altitude, the change in genetic distance
# per unit (here one SD) of altitude difference, conditional on resistance.
summary(fit_joint)
#> Conductance surface with 2531 vertices (31 focal) estimated by maximum likelihood
#> Call:   terradish(formula = Fst ~ altitude + forestcover, data = surface, 
#>     measurement_model = g_joint)
#> 
#> Loglikelihood: 1541.802 (7 degrees freedom)
#> AIC: -3069.603 
#> 
#> Number of function calls: 18 
#> Number of optimization steps: 5 
#> Convergence: converged | projected gradient: 4.2e-09 
#> Norm of gradient at MLE: 5.801505e-09 
#> 
#> Nuisance parameters (conditional on fitted conductance surface; measurement-model scale):
#>                                  Estimate   Std. Error  Lower 95%  Upper 95%
#> alpha                            -0.055479   0.017745   -0.090258  -0.020699
#> beta                              0.387706   0.013284    0.361669   0.413743
#> absdiff_spatial_template_lowres   0.014297   0.001392    0.011569   0.017024
#> tau                               5.323936   0.243868    4.845964   5.801908
#> rho                               2.587655   0.271940    2.054663   3.120647
#> 
#> Resistance-distance equivalent of one unit of environmental difference:
#>                        covariate      ratio          SE      lower      upper
#>  absdiff_spatial_template_lowres 0.03687554 0.005040606 0.02699613 0.04675495
#>             family
#>  gaussian_distance
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> altitude    -0.83600    0.10940  -7.642 2.15e-14 ***
#> forestcover -0.07662    0.17215  -0.445    0.656    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Correlation of Coefficients:
#>              altitude
#> forestcover 0.7029914
terradish_ibe_ratio(fit_joint)
#>                         covariate      ratio          SE      lower      upper
#> 1 absdiff_spatial_template_lowres 0.03687554 0.005040606 0.02699613 0.04675495
#>              family
#> 1 gaussian_distance
# }
```
