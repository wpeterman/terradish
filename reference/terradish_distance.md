# Resistance distances on a grid of conductance parameters

Calculates resistance distances associated with a parameterized
conductance surface across a grid of parameter values.

## Usage

``` r
terradish_distance(
  theta,
  formula,
  data,
  conductance_model = loglinear_conductance,
  conductance = TRUE,
  covariance = FALSE,
  cores = 1L
)
```

## Arguments

- theta:

  A matrix of dimension (grid size) x (number of parameters)

- formula:

  A formula whose right-hand side names the covariates used to create
  `data`. A left-hand side, if supplied, is not used, so a one-sided
  formula such as `~ forestcover + altitude` is sufficient.

- data:

  An object of class `terradish_graph` (see
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md))

- conductance_model:

  A function of class `terradish_conductance_model_factory` (see
  [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md))

- conductance:

  Retained for backward compatibility. Only `conductance = TRUE` is
  currently implemented.

- covariance:

  If `TRUE`, instead of a matrix of resistance distances, return the
  associated submatrix of the generalized inverse of graph Laplacian

- cores:

  Number of worker processes to use. `1` evaluates the grid serially.

## Value

A list of class `c("terradish_grid", "radish_grid")`:

- theta:

  Data frame of the evaluated conductance coefficients, one row per grid
  point.

- distance:

  When `covariance = FALSE`, a focal-site by focal-site by grid-point
  array of resistance distances. Larger values mean greater graph
  resistance between two sites. Their units depend on the conductance
  scale, so compare them within a grid rather than as map distances.

- covariance:

  When `covariance = TRUE`, the corresponding array of graph covariance
  (the generalized inverse of the graph Laplacian), returned in place of
  `distance`.

## See also

[`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other model fitting:
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md),
[`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_algorithm()`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md),
[`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)

## Examples

``` r
library(terra)

data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)[1:10]
names(covariates) <- c("altitude", "forestcover")
# Coarsened so the example runs quickly; name layers, then scale.
covariates <- scale_covariates(terra::aggregate(covariates, fact = 5,
                                                na.rm = TRUE))
# Keep one site per raster cell.
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
surface <- conductance_surface(covariates, coords[keep])

theta <- as.matrix(expand.grid(forestcover = seq(-0.5, 0.5, length.out = 3),
                               altitude = seq(-0.5, 0.5, length.out = 3)))

distances <- terradish_distance(theta, ~ forestcover + altitude, surface)

# Resistance distances under uniform conductance (all coefficients zero)
ibd <- which(theta[, 1] == 0 & theta[, 2] == 0)
round(distances$distance[, , ibd], 4)
#>         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]  [,10]
#>  [1,] 0.0000 0.1216 0.3890 0.3269 0.3312 0.1656 0.2103 0.1589 0.4220 0.4402
#>  [2,] 0.1216 0.0000 0.3796 0.3179 0.3220 0.1875 0.1955 0.1840 0.4127 0.4310
#>  [3,] 0.3890 0.3796 0.0000 0.2512 0.2363 0.4136 0.3465 0.4182 0.2233 0.2642
#>  [4,] 0.3269 0.3179 0.2512 0.0000 0.1535 0.3500 0.2867 0.3556 0.2706 0.2827
#>  [5,] 0.3312 0.3220 0.2363 0.1535 0.0000 0.3550 0.2897 0.3602 0.2629 0.2793
#>  [6,] 0.1656 0.1875 0.4136 0.3500 0.3550 0.0000 0.2437 0.1253 0.4462 0.4642
#>  [7,] 0.2103 0.1955 0.3465 0.2867 0.2897 0.2437 0.0000 0.2449 0.3802 0.3989
#>  [8,] 0.1589 0.1840 0.4182 0.3556 0.3602 0.1253 0.2449 0.0000 0.4510 0.4692
#>  [9,] 0.4220 0.4127 0.2233 0.2706 0.2629 0.4462 0.3802 0.4510 0.0000 0.1890
#> [10,] 0.4402 0.4310 0.2642 0.2827 0.2793 0.4642 0.3989 0.4692 0.1890 0.0000
```
