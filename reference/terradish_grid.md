# Profile likelihood on a grid of conductance parameters

Calculates the profile likelihood of a parameterized conductance surface
across a grid of parameter values (that is, the nuisance parameters are
optimized at each point on the grid).

## Usage

``` r
terradish_grid(
  theta,
  formula,
  data,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  nu = NULL,
  nonnegative = TRUE,
  conductance = TRUE,
  covariance = FALSE,
  cores = 1L,
  approximation = "none",
  approximation_control = NULL
)
```

## Arguments

- theta:

  A matrix of dimension (grid size) x (number of parameters)

- formula:

  A formula with the name of a response matrix compatible with the
  selected measurement model on the left-hand side, and covariates used
  to create `data` on the right-hand side.

- data:

  An object of class `terradish_graph` (see
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md))

- conductance_model:

  A function of class `terradish_conductance_model_factory` (see
  [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md))

- measurement_model:

  A function of class `terradish_measurement_model` (see
  [`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md))

- nu:

  Effective Wishart degrees of freedom passed to measurement models that
  require it. This is a user-supplied effective-information parameter,
  not the marker count. Report sensitivity to plausible values; see
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).

- nonnegative:

  Force regression-like `measurement_model` to have nonnegative slope?

- conductance:

  Retained for backward compatibility. Only `conductance = TRUE` is
  currently implemented.

- covariance:

  If `TRUE`, additionally return (a submatrix of) the generalized
  inverse of graph Laplacian across the grid

- cores:

  Number of worker processes to use. `1` evaluates the grid serially.

- approximation:

  Only `"none"` is supported; every grid point uses the full graph and
  focal set.

- approximation_control:

  Reserved for compatibility; unused.

## Value

A list of class `c("terradish_grid", "radish_grid")`, with one entry per
row of `theta`:

- theta:

  Data frame of the evaluated conductance coefficients (log conductance
  per unit of each covariate).

- loglik:

  Numeric vector of profile log-likelihoods, with nuisance parameters
  optimized at each grid point. Larger values indicate better fit, and
  differences between grid points are on the log-likelihood scale.

- phi:

  Matrix of optimized nuisance parameters, one column per grid point, on
  the measurement model's internal scale.

- covariance:

  When `covariance = TRUE`, a focal-site by focal-site by grid-point
  array of graph covariance (the generalized inverse of the graph
  Laplacian); otherwise `NULL`.

- approximation:

  Record of the approximation setting (always `"none"`).

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
for the optimizer and
[`terradish_distance`](https://wpeterman.github.io/terradish/reference/terradish_distance.md)
for resistance distances on the same kind of grid.

Other model fitting:
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md),
[`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md),
[`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md),
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)

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
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords[keep])

fit_mlpe <- terradish(Fst ~ altitude + forestcover, data = surface,
                      measurement_model = mlpe)

theta <- as.matrix(expand.grid(forestcover = seq(-0.5, 0.5, length.out = 3),
                               altitude = seq(-0.5, 0.5, length.out = 3)))

grid <- terradish_grid(theta, Fst ~ forestcover + altitude, surface,
                       measurement_model = mlpe)

# Profile log-likelihood at each grid point; the maximum should lie near
# the fitted coefficients.
cbind(grid$theta, loglik = grid$loglik)
#>   forestcover altitude   loglik
#> 1        -0.5     -0.5 1468.833
#> 2         0.0     -0.5 1492.091
#> 3         0.5     -0.5 1485.862
#> 4        -0.5      0.0 1403.162
#> 5         0.0      0.0 1443.254
#> 6         0.5      0.0 1460.622
#> 7        -0.5      0.5 1344.748
#> 8         0.0      0.5 1384.486
#> 9         0.5      0.5 1401.381
coef(fit_mlpe)
#>    altitude forestcover 
#> -0.79206795 -0.09612631 
# }
```
