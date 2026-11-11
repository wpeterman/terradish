# Change in genetic distance across an MLPE pairwise covariate contrast

Computes the expected change in the genetic-distance response associated
with changing one or more pairwise covariates in an
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
measurement model. This is useful for interpreting IBE coefficients,
which are profiled nuisance parameters rather than conductance-surface
coefficients.

## Usage

``` r
mlpe_response_change(
  object,
  covariate = NULL,
  probs = c(0.1, 0.9),
  values = NULL,
  conf.level = 0.95
)
```

## Arguments

- object:

  A fitted `terradish` model whose measurement model was created with
  [`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md).

- covariate:

  Optional character vector naming the pairwise covariates to summarize.
  Defaults to all pairwise covariates stored in the measurement model.

- probs:

  Numeric vector of length two giving the lower and upper quantiles used
  to define the contrast when `values` is `NULL`. Defaults to the 10th
  and 90th percentiles.

- values:

  Optional numeric vector of length two giving explicit lower and upper
  covariate values for a single `covariate`. Alternatively, provide a
  named list with one length-two numeric vector per covariate.

- conf.level:

  Confidence level for the Wald confidence interval.

## Value

A data frame with one row per covariate and columns describing the low
and high covariate values, the contrast, estimated response-scale
change, confidence limits, and confidence level.

## Details

For a model fitted with
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
the pairwise mean structure includes terms of the form \\\gamma
Z\_{ij}\\. This function reports \\\gamma (z\_{\mathrm{high}} -
z\_{\mathrm{low}})\\, the expected response-scale change in genetic
distance across the requested covariate contrast, conditional on the
fitted conductance surface.

The interval uses the conditional standard error for the MLPE covariate
coefficient from `summary(object)` and treats the selected covariate
contrast as fixed. The reported contrast is a conditional association
under the fitted MLPE model. It is not a causal effect or proof that IBE
has been separated from IBR or other spatially correlated predictors.

## Worked examples

[`vignette("v08-ibe-ibr-workflow", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v08-ibe-ibr-workflow.md).

## See also

[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_endpoint_covariates`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

## Examples

``` r
# \donttest{
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

# site-level altitude, turned into a pairwise |difference| IBE predictor
site_altitude <- data.frame(
  altitude = terra::extract(terra::unwrap(melip.altitude), coords)[, 2])
g_joint <- mlpe_covariates(site_altitude)

# joint isolation by environment and isolation by resistance
fit_joint <- terradish(Fst ~ altitude + forestcover, data = surface,
                       conductance_model = loglinear_conductance,
                       measurement_model = g_joint)

# Expected change in genetic distance across the covariate's 10th-to-90th
# percentile spread, on the response scale rather than the log scale.
change <- mlpe_response_change(fit_joint)
change
#>          covariate        low      high  contrast   estimate   conf.low
#> 1 absdiff_altitude 0.02037562 0.5275703 0.5071947 0.03288384 0.02661041
#>    conf.high conf.level
#> 1 0.03915727       0.95

# A narrower, interquartile contrast
mlpe_response_change(fit_joint, probs = c(0.25, 0.75))
#>          covariate     low      high  contrast   estimate   conf.low  conf.high
#> 1 absdiff_altitude 0.08377 0.3477135 0.2639435 0.01711271 0.01384802 0.02037739
#>   conf.level
#> 1       0.95
# }
```
