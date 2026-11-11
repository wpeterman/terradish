# Least-squares measurement model

A function of class `"terradish_measurement_model"` that evaluates a
Gaussian likelihood with independent errors and a linear mean structure
relating observed genetic distances to resistance distances. This is the
fastest measurement model and a useful exploratory or diagnostic choice.
For distance-response inference,
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md)
represents shared-site dependence that this model omits.

## Usage

``` r
leastsquares(
  E,
  S,
  phi,
  nu = NULL,
  gradient = TRUE,
  hessian = TRUE,
  partial = TRUE,
  nonnegative = TRUE,
  validate = FALSE
)
```

## Arguments

- E:

  Conductance-implied covariance matrix: the generalized inverse of the
  graph Laplacian, evaluated at the current conductance parameters.
  Passed automatically by the optimizer; users normally do not call this
  function directly.

- S:

  Square, symmetric matrix of observed pairwise genetic distances (e.g.
  \\F\_{ST}\\). Must have the same dimensions as `E`.

- phi:

  Named numeric vector of nuisance parameters `(alpha, beta, tau)`. Omit
  to obtain maximum-likelihood starting values from an
  [`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html) fit.

- nu:

  Unused; present for a common interface with Wishart measurement
  models.

- gradient:

  Logical. Compute gradient of the negative log-likelihood with respect
  to `phi`?

- hessian:

  Logical. Compute Hessian of the negative log-likelihood with respect
  to `phi`?

- partial:

  Logical. Compute second partial derivatives with respect to `phi`,
  `E`, and `S`? Required by the optimizer; set `FALSE` only for
  standalone likelihood evaluations.

- nonnegative:

  Logical. Constrain the IBR slope `beta` to be nonnegative? Default
  `TRUE` prevents a nonsensical negative resistance-distance slope.

- validate:

  Logical. Numerically validate gradients and Hessians via numDeriv?
  Very slow; intended for debugging small examples only.

## Value

When `phi` is missing, a list with elements `phi` (starting values),
`lower`, and `upper` (parameter bounds). Otherwise a list containing:

- objective:

  Negative log-likelihood.

- fitted:

  Matrix of expected genetic distances (same dimensions as `S`).

- boundary:

  Logical; `TRUE` if the MLE is on the boundary (`beta = 0`), indicating
  no detectable IBR signal.

- gradient:

  Gradient of the negative log-likelihood with respect to `phi` (if
  `gradient = TRUE`).

- hessian:

  Hessian matrix with respect to `phi` (if `hessian = TRUE`).

- gradient_E:

  Gradient with respect to `E` (if `partial = TRUE`).

- partial_E:

  Jacobian of `gradient_E` with respect to `phi` (if `partial = TRUE`).

- partial_S:

  Jacobian of `gradient` with respect to the lower triangle of `S` (if
  `partial = TRUE`).

- jacobian_E:

  Function for reverse-mode algorithmic differentiation through `E` (if
  `partial = TRUE`).

- jacobian_S:

  Function for reverse-mode algorithmic differentiation through `S` (if
  `partial = TRUE`).

## Details

The nuisance parameters are:

- `alpha`:

  Intercept of the regression of genetic distance on resistance
  distance.

- `beta`:

  IBR slope; constrained \\\geq 0\\ when `nonnegative = TRUE`.

- `tau`:

  Log-precision of the Gaussian errors: residual variance is
  \\\exp(-\tau)\\, so the residual standard deviation is
  `exp(-tau / 2)`. In the Wishart measurement models `tau` denotes a
  different quantity (see
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)).

The mean structure is \\S\_{ij} = \alpha + \beta R\_{ij} + e\_{ij}\\,
where \\R\_{ij}\\ is the resistance distance between sites \\i\\ and
\\j\\ (derived from `E`) and \\e\_{ij}\\ are independent Gaussian errors
with common precision \\\exp(\tau)\\.

Pairwise genetic distances are not independent: any two pairs sharing a
sampling site can be correlated. `leastsquares` ignores this dependence,
so its standard errors can be too small. Use it for fast exploration,
diagnostics, or initial parameter searches rather than as the default
basis for distance-response inference.

## See also

[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)

## Examples

``` r
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

# Graph covariance (generalized inverse of the graph Laplacian) among the
# focal sites under uniform conductance (all coefficients zero).
laplacian_inv <- terradish_distance(theta = matrix(0, 1, 2),
                                    formula = ~ forestcover + altitude,
                                    data = surface,
                                    covariance = TRUE)$covariance[, , 1]

# Without 'phi': starting values and bounds for the nuisance parameters
leastsquares(laplacian_inv, Fst)
#> $phi
#>      alpha       beta        tau 
#> 0.01002299 0.22923592 6.47733469 
#> 
#> $lower
#> [1] -Inf    0 -Inf
#> 
#> $upper
#> [1] Inf Inf Inf
#> 
# With 'phi': the negative log-likelihood and derivatives at those values
value <- leastsquares(laplacian_inv, Fst,
                      phi = c(alpha = 0, beta = 0.5, tau = -0.1))
value$objective
#> [1] 25.83699
```
