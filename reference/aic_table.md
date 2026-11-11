# Rank fitted models by AIC, AICc, or BIC

Creates a model-selection table from fitted `terradish` models using
AIC, AICc, or BIC.

## Usage

``` r
aic_table(
  mod_list,
  AICc = FALSE,
  BIC = FALSE,
  mod_names = NULL,
  verbose = FALSE
)
```

## Arguments

- mod_list:

  List of fitted `terradish` models.

- AICc:

  Should second-order AIC (Akaike's information criterion with
  small-sample correction) be used instead of AIC? When `TRUE`, the
  correction uses \\n\\ = number of focal sampling sites (not the number
  of pairwise observations). A common guideline is to prefer AICc over
  AIC when \\n / K \< 40\\, where \\K\\ is the number of estimated
  parameters.

- BIC:

  Should BIC be used instead of AIC? When `TRUE`, the penalty uses \\n\\
  = number of focal sampling sites (individuals or populations), the
  same sample size as AICc, not the number of pairwise observations.
  Pairs that share a site are not independent, so the number of pairs
  overstates the information in the data.

- mod_names:

  Optional model names. By default the right-hand side of each fitted
  formula is used. For MLPE measurement models with additional pairwise
  covariates, the default appends `[mlpe:n]`, where `n` is the number of
  added pairwise covariate columns.

- verbose:

  Should the table be printed to the console?

## Value

A data frame with one row per model, sorted from best to worst: `model`
(the label), `K` (number of estimated parameters, conductance plus
nuisance), the criterion (`AIC`, `AICc`, or `BIC`; smaller is better),
its difference from the best model (`Delta_AIC` and so on), the model
weight and cumulative weight, and `loglik`. The `converged` column is
true only when both the outer optimizer and final nuisance optimization
record code 0. Missing convergence diagnostics in older objects are
treated as unverified (false); `boundary` marks active parameter bounds
or no resistance structure. A warning identifies tables containing
nonconverged or boundary fits.

## Details

At nominal Wishart information, information criteria can favor
extensions that are unsupported across replicate histories. Use
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
to select conductance terms and
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
to assess inferential sensitivity.

Information-criterion comparison requires the same observed response,
focal sites, graph domain, and likelihood family. Wishart fits must also
use the same effective degrees of freedom, `nu`. `leastsquares` and
`mlpe` are Gaussian likelihoods for the same pairwise-distance response
and can be ranked when those conditions hold. `generalized_wishart` and
`wishart_covariance` use different response representations and must not
be ranked together. Cross-family rankings, such as `mlpe` versus
`generalized_wishart`, are invalid.

The function checks response values, model family, fitted dimensions,
recorded `nu`, and retained graph signatures. Older fits without
signatures require the user to verify the graph domain.

## Worked example

See
[`vignette("v04-model-comparison", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v04-model-comparison.md)
for compatible model rankings, coefficient plots, and predictive
comparisons.

## See also

[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
for predictive comparison,
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
for sensitivity to `nu`, and
[`anova()`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)
for nested likelihood-ratio tests.

Other inference and model comparison:
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_folds()`](https://wpeterman.github.io/terradish/reference/terradish_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

## Examples

``` r
# \donttest{
library(terra)
#> terra 1.9.34

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

fit1 <- terradish(Fst ~ altitude + forestcover, surface,
                  measurement_model = mlpe)
fit2 <- terradish(Fst ~ altitude * forestcover, surface,
                  measurement_model = mlpe)

# One row per model, best first. K counts conductance plus nuisance
# parameters; Delta_AIC is the difference from the best model.
aic_table(list(fit1, fit2))
#>                                           model K       AIC Delta_AIC AIC_wt
#> 2 altitude + forestcover + altitude:forestcover 7 -3007.215    0.0000      1
#> 1                        altitude + forestcover 6 -2978.092   29.1238      0
#>   Cum.wt   loglik converged boundary
#> 2      1 1510.608      TRUE    FALSE
#> 1      1 1495.046      TRUE    FALSE
aic_table(list(fit1, fit2), AICc = TRUE)
#>                                           model K       AIC      AICc
#> 2 altitude + forestcover + altitude:forestcover 7 -3007.215 -3002.346
#> 1                        altitude + forestcover 6 -2978.092 -2974.592
#>   Delta_AICc AICc_wt Cum.wt   loglik converged boundary
#> 2     0.0000       1      1 1510.608      TRUE    FALSE
#> 1    27.7542       0      1 1495.046      TRUE    FALSE
# }
```
