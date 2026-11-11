# Maximum likelihood population effects (MLPE) measurement model

A function of class `"terradish_measurement_model"` that evaluates the
maximum likelihood population effects (MLPE) likelihood. MLPE accounts
for the non-independence of pairwise genetic distances that share a
sampling site. It is the preferred distance-response likelihood when a
defensible effective Wishart degrees-of-freedom value is unavailable or
the response is not admissible for
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md).

## Usage

``` r
mlpe(
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

  Named numeric vector of nuisance parameters `(alpha, beta, tau, rho)`.
  The stored `tau` is log precision and the stored `rho` is an
  unconstrained, logit-scale parameter. Omit to obtain starting values
  from a least-squares fit and a weak positive MLPE correlation.

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
  `E`, and `S`? Required by the optimizer.

- nonnegative:

  Logical. Constrain the IBR slope `beta` to be nonnegative? Default
  `TRUE`.

- validate:

  Logical. Numerically validate gradients and Hessians via numDeriv?
  Very slow; for debugging small examples only.

## Value

When `phi` is missing, a list with elements `phi` (starting values),
`lower`, and `upper` (parameter bounds). Otherwise a list containing:

- objective:

  Negative log-likelihood.

- fitted:

  Matrix of expected genetic distances (same dimensions as `S`).

- boundary:

  Logical; `TRUE` if `beta = 0`, indicating no detectable IBR signal.

- gradient:

  Gradient with respect to `phi` (if `gradient = TRUE`).

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

  Function for reverse-mode AD through `E` (if `partial = TRUE`).

- jacobian_S:

  Function for reverse-mode AD through `S` (if `partial = TRUE`).

## Details

The likelihood uses maximum likelihood (ML), not restricted maximum
likelihood (REML). AIC can therefore compare different conductance
formulas fitted to the same response, sites, graph, and likelihood
family. That comparison does not validate the correlation model: in
unpublished development simulations, model-based standard errors
understated replicate variation about 1.7-fold. Use spatial CV to assess
out-of-sample support for conductance terms.

With the standard nuisance bounds, fitting profiles the regression
coefficients and residual precision analytically at each shared-site
correlation, then searches the correlation before Newton refinement. The
resistance predictor is centered and scaled internally so changing its
units does not change the nuisance optimum or convergence criterion.
Stored coefficients remain in the original units. Fixed predictive
parameters and custom nuisance bounds use their supplied constraints. A
profile reaching the upper numerical correlation endpoint is rejected
because a regular finite nuisance optimum has not been established.

[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
reports `rho` on its internal logit scale, with a Wald interval on that
scale. Convert an estimate or interval endpoint with `plogis(rho) / 2`
to obtain the shared-site correlation, which lies between 0 and 0.5. A
value near 0.5 means site identity explains most of the residual
variation. When the correlation approaches zero (fitted internal logit
below -8), the logit parameter is poorly identified, so
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
reports `rho` as 0 with no standard error and prints a boundary note.
The parameter remains in the degrees of freedom because it was
estimated; fixing it in advance would define a different model. This
boundary differs from a zero resistance coefficient, which leaves
conductance unidentified.

To read the other nuisance parameters, `alpha` is the expected genetic
distance at zero resistance distance, `beta` is the increase in genetic
distance per unit of resistance distance, and the residual standard
deviation is `exp(-tau / 2)`.

The nuisance parameters are:

- `alpha`:

  Intercept of the mean structure.

- `beta`:

  Resistance-distance (IBR) slope; constrained \\\geq 0\\ when
  `nonnegative = TRUE`.

- `tau`:

  Log-precision parameter: residual variance is \\\exp(-\tau)\\. In the
  Wishart measurement models the same name denotes a different quantity,
  the multiplier on graph covariance (see
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)).

- `rho`:

  Unconstrained optimization-scale parameter. The actual shared-site
  correlation is \\\rho\_{MLPE} = \text{plogis}(\rho)/2 \in (0, 0.5)\\.

The mean structure is \\S\_{ij} = \alpha + \beta R\_{ij} + e\_{ij}\\,
where \\R\_{ij}\\ is the resistance distance derived from `E`. The
residual vector \\e\\ follows the MLPE correlation structure of Clarke
et al. (2002): two pairs \\(i,j)\\ and \\(i,k)\\ sharing site \\i\\ have
correlation \\\rho\_{MLPE}\\, while pairs sharing no site are
uncorrelated.

Use
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md)
to extend the mean structure with additional fixed pairwise covariates
(isolation by environment). Coefficients describe conditional
associations under the specified mean and correlation structure; they do
not isolate causal environmental effects.

## Worked example

[`vignette("v03-mlpe-pairwise-distance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v03-mlpe-pairwise-distance.md)
fits and interprets MLPE models for a pairwise \\F\_{ST}\\ matrix.

## References

Clarke RT, Rothery P, Raybould AF (2002). Confidence limits for
regression relationships between distance matrices: estimating gene flow
with distance. Journal of Agricultural, Biological, and Environmental
Statistics, 7, 361–372.
[doi:10.1198/108571102320](https://doi.org/10.1198/108571102320) .

## See also

[`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
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
mlpe(laplacian_inv, Fst)
#> $phi
#>       alpha        beta         tau         rho 
#>  0.01002299  0.22923592  6.47733469 -1.38629436 
#> 
#> $lower
#> [1] -Inf    0 -Inf -Inf
#> 
#> $upper
#> [1] Inf Inf Inf Inf
#> 
# With 'phi': the negative log-likelihood and derivatives at those values.
# qlogis(2 * 0.2) sets the actual shared-site correlation to 0.2.
value <- mlpe(laplacian_inv, Fst,
              phi = c(alpha = 0, beta = 0.5, tau = -0.1,
                      rho = qlogis(2 * 0.2)))
value$objective
#>       tau 
#> -58.01211 
```
