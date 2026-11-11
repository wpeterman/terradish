# Search raster smoothing scales by refitting the model

Performs an outer optimization over raster scale parameters by
repeatedly rescaling raster covariates, rebuilding or updating the
conductance surface, and refitting
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Usage

``` r
terradish_scale_optim(
  formula,
  covariates,
  coords,
  scales = NULL,
  lower = NULL,
  upper = NULL,
  scale_grid = NULL,
  search = c("coordinate", "grid"),
  objective = c("aic", "logLik"),
  grid_points = 5L,
  maxit = 10L,
  tol = 0.1,
  scale_fun = "terradish",
  scale_args = NULL,
  postprocess = scale_covariates,
  directions = 8L,
  theta = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- formula:

  A terradish model formula.

- covariates:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  containing the raster layers whose smoothing scales will be optimized.

- coords:

  Focal coordinates passed to
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- scales:

  Optional starting values for the scale parameters. Required for
  `search = "coordinate"` unless `scale_grid` is supplied.

- lower, upper:

  Lower and upper bounds for each scale parameter, in map units for both
  built-in smoothing methods. The default uses the same normalized,
  truncated Gaussian as the joint scale model. A useful upper bound is
  one sixth of the smaller raster dimension times its cell width. Used
  to construct `scale_grid` when that is not supplied, and required for
  `search = "coordinate"`.

- scale_grid:

  Optional named list of candidate scale values for `search = "grid"`.

- search:

  One of `"coordinate"` or `"grid"`.

- objective:

  Objective used for the outer optimization: `"aic"` or `"logLik"`.

- grid_points:

  Number of equally spaced candidate scales per raster when
  `search = "grid"` and `scale_grid` is not supplied.

- maxit:

  Maximum number of outer coordinate-search iterations.

- tol:

  Convergence tolerance for coordinate search, measured as the maximum
  absolute change in the scale vector between iterations.

- scale_fun:

  Smoothing method: `"terradish"` (default) or `"multiScaleR"`. Both use
  map units, but their kernel constructions differ, so their fitted
  scales are not interchangeable. A custom function may instead accept a
  single-layer raster and a numeric `scale=` argument. `NULL` selects
  the default.

- scale_args:

  Optional named list passed to `scale_fun`.

- postprocess:

  Optional function applied to the full scaled raster stack before
  constructing the conductance surface. Defaults to
  [`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)
  so each scaled layer is re-standardized. Set to `NULL` to use the
  scaled rasters as returned by `scale_fun`.

- directions:

  Neighborhood definition passed to
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- theta:

  Optional starting values passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
  at each outer evaluation.

- verbose:

  Print outer-optimization progress?

- ...:

  Additional arguments passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Value

An object of class `"terradish_scale_optim"`, a list with:

- par:

  The selected scale for each raster layer, in the units `scale_fun`
  uses (map units for both built-in methods).

- value:

  The outer objective at `par`. Lower is better, whichever objective was
  used.

- fit:

  The fitted `terradish` model at `par`. Read its conductance
  coefficients with [`coef()`](https://rdrr.io/r/stats/coef.html) the
  same way as any other fit; they are on the log-conductance scale.

- surface:

  The conductance surface built from the rasters smoothed at `par`,
  ready to pass to
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  or [`plot()`](https://rspatial.github.io/terra/reference/plot.html).

- evaluations:

  A data frame with one row per candidate scale combination: the scale
  for each layer, the `objective` value, `aic`, `loglik`, and a `status`
  column that records `"OK"` or the error for candidates that failed to
  fit. Inspect this to see how sharply the objective distinguishes the
  scales; a flat column means the objective does not identify a
  smoothing scale within the searched range.

- search:

  The search strategy used, `"grid"` or `"coordinate"`.

- objective:

  The name of the objective minimized, `"aic"` or `"logLik"`.

- iterations, convergence:

  Number of outer iterations, and `0` if the coordinate search converged
  (always `0` for a grid search).

## Details

This helper treats scale optimization as an outer problem around the
usual
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
fit, rather than as an inner Newton/BFGS parameter. That makes it
compatible with raster-scaling approaches such as
[`multiScaleR::kernel_scale.raster()`](https://rdrr.io/pkg/multiScaleR/man/kernel_scale.raster.html),
while keeping the terradish core optimizer unchanged. The selected
values describe the transformation performed by `scale_fun` in that
function's units. They are tuning parameters for the raster
representation, not direct estimates of dispersal distance, home-range
size, organismal movement, or a causal ecological-process scale. Report
the search range and objective, and account for this outer selection
when describing model-selection uncertainty. Each searched scale counts
toward the parameter total used for information criteria. Prefer the
joint Gaussian model when its kernel matches the scientific question;
this outer search is a secondary option for other kernels or an explicit
scale grid. The default terradish kernel matches the joint model,
whereas multiScaleR uses a different construction even when both scales
are expressed in map units.

Coordinate search updates one raster scale at a time using
[`optimize`](https://rdrr.io/r/stats/optimize.html) while holding the
others fixed. Grid search evaluates every combination in `scale_grid`.
In both modes, failed terradish fits receive an infinite outer objective
and are retained in the returned evaluation table.

## See also

[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md)
for joint estimation of Gaussian scales,
[`gaussian_scale_summary`](https://wpeterman.github.io/terradish/reference/gaussian_scale_summary.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other spatial scale:
[`gaussian_scale_profile()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md),
[`gaussian_scale_summary()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_summary.md),
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md)

## Examples

``` r
# \donttest{
data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)
# Name the layers. Do not scale them here: the default `postprocess`
# re-standardizes each layer after it is smoothed.
names(covariates) <- c("altitude", "forestcover")

# Every candidate scale re-smooths the rasters and refits, so coarsen the
# grid for the example. Use the full resolution in a real analysis.
covariates <- terra::aggregate(covariates, fact = 3, na.rm = TRUE)

# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
Fst <- melip.Fst[keep, keep]

# Express the bounds in map units, using the coarsened cell width.
cell_width <- mean(terra::res(covariates))
# A 2-point grid over two layers is 4
# candidate fits; use a wider grid in practice.
fit_scale <- terradish_scale_optim(
  Fst ~ altitude + forestcover,
  covariates = covariates,
  coords = coords[keep],
  lower = cell_width,
  upper = 6 * cell_width,
  search = "grid",
  grid_points = 2,
  measurement_model = leastsquares,
  optimizer = "bfgs"
)
fit_scale$par          # selected smoothing scale per layer, in map units
#>    altitude forestcover 
#>       0.135       0.135 
fit_scale$value        # objective value there (lower is better)
#> [1] -2748.498
fit_scale$objective    # which objective was minimized ("aic" by default)
#> [1] "aic"
fit_scale$evaluations  # every candidate tried, with its AIC and loglik
#>   altitude forestcover objective       aic   loglik status
#> 1    0.135       0.135 -2748.498 -2748.498 1381.249     OK
#> 2    0.810       0.135 -2743.431 -2743.431 1378.716     OK
#> 3    0.135       0.810 -2728.719 -2728.719 1371.360     OK
#> 4    0.810       0.810 -2726.710 -2726.710 1370.355     OK
coef(fit_scale$fit)    # conductance coefficients at the selected scale
#>    altitude forestcover 
#>  -0.2419000   0.2382399 
# }
```
