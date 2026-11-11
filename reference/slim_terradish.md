# Shrink a fitted model for storage

Removes components that are useful for prediction or diagnostics but can
make a fitted `terradish` object expensive to serialize. The fitted
coefficients, likelihood, Hessian, fitted response, convergence
diagnostics, and approximation metadata are retained unless explicitly
dropped.

## Usage

``` r
slim_terradish(object, drop = c("submodels", "leverage"))
```

## Arguments

- object:

  A fitted `terradish` object.

- drop:

  Character vector naming optional components to remove. Supported
  values are `"submodels"`, `"response"`, and `"leverage"`. Dropping
  `"submodels"` removes the fitted conductance and measurement-model
  closures, which usually account for most of the serialized object
  size. Dropping `"response"` also removes the response copy stored in
  the comparison contract.

## Value

A fitted `terradish` object with a `storage` element that records the
removed components and the object's size before and after slimming.

## Details

A slim object remains suitable for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rspatial.github.io/terra/reference/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html), and fitted-value extraction
and residual-permutation simulation. Operations that must reevaluate the
conductance or measurement model, including conductance- surface
prediction and most plot types, require the retained model closures and
will stop with an informative error if those closures were dropped.

Slimming is irreversible. Retain the original fit, model formula, graph,
and response when later prediction, simulation, or refitting may be
needed.

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
(argument `slim`),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md)

Other solver settings:
[`print.terradish_setting_assessment()`](https://wpeterman.github.io/terradish/reference/print.terradish_setting_assessment.md),
[`terradish_assess_settings()`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md),
[`terradish_solver_benchmark()`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)

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
fit <- terradish(Fst ~ altitude + forestcover, surface,
                 measurement_model = mlpe)
fit_slim <- slim_terradish(fit)
fit_slim$storage[c("bytes_before", "bytes_after")]
#> $bytes_before
#> [1] 1508584
#> 
#> $bytes_after
#> [1] 49512
#> 
coef(fit_slim)     # estimates, summaries, and intervals still work
#>    altitude forestcover 
#> -0.79206795 -0.09612631 
confint(fit_slim)
#>                  2.5 %     97.5 %
#> altitude    -1.0862013 -0.4979346
#> forestcover -0.4735909  0.2813383
# }
```
