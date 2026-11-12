# Generalized Wishart measurement model for squared distances

A function of class `"terradish_measurement_model"` that evaluates the
generalized Wishart likelihood for an admissible pairwise genetic
squared-distance matrix. This model requires a distance representation
that is coherently related to a centered positive semidefinite
covariance matrix; it is not valid for an arbitrary dissimilarity
matrix. The effective marker degrees of freedom must also be supplied.
For covariance-matrix data, use
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).

## Usage

``` r
generalized_wishart(
  E,
  S,
  phi,
  nu,
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
  graph Laplacian at the current conductance parameters. Passed
  automatically by the optimizer.

- S:

  Square, symmetric admissible genetic squared-distance matrix with zero
  diagonal. It must have the same dimensions as `E` and must correspond
  to a centered positive semidefinite covariance representation under
  the scaling used to construct it.

- phi:

  Named numeric vector of nuisance parameters `(tau, sigma)`. Omit to
  obtain default starting values `c(1, 0)`.

- nu:

  One finite positive number describing effective Wishart information.
  You must supply it; the fit does not estimate it. It is not the marker
  count. See the effective-information section and report sensitivity
  across plausible values with
  [`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md).

- gradient:

  Logical. Compute gradient of the negative log-likelihood with respect
  to `phi`?

- hessian:

  Logical. Compute Hessian with respect to `phi`?

- partial:

  Logical. Compute second partial derivatives with respect to `phi`,
  `E`, and `S`?

- nonnegative:

  Unused; present for interface consistency.

- validate:

  Logical. Numerically validate gradients and Hessians via numDeriv?
  Very slow; for debugging small examples only.

## Value

When `phi` is missing, a list with elements `phi` (named starting values
`c(tau = 1, sigma = 0)`), `lower` (`c(0, -Inf)`), and `upper`. Otherwise
a list containing:

- objective:

  Negative log-likelihood.

- fitted:

  Matrix of expected squared genetic distances.

- boundary:

  Logical; `TRUE` if `tau = 0`, indicating no detectable IBR signal.

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

The nuisance parameters are:

- `tau`:

  Nonnegative scale applied to the conductance-implied covariance `E`. A
  `tau` near zero signals no detectable contribution from
  conductance-implied covariance under this model. This differs from
  `tau` in
  [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) and
  [`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
  where it is the log residual precision.

- `sigma`:

  Log-scale identity component: the nugget variance added to the
  diagonal is \\\exp(\sigma)\\. It represents diagonal variance not
  captured by the graph covariance, including sampling variation and
  model mismatch. It is not a smoothing scale or a direct estimate of
  genetic drift or population size. The simulators
  [`simulate_covariance_response`](https://wpeterman.github.io/terradish/reference/simulate_covariance_response.md)
  and
  [`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md)
  take the nugget on its natural scale, so pass them
  [`exp()`](https://rdrr.io/r/base/Log.html) of a fitted `sigma`.

The fitted covariance is \\\Sigma = \tau E + \exp(\sigma) I\\. The model
evaluates the generalized Wishart log-likelihood for the observed
admissible squared-distance matrix `S` after projecting out the grand
mean, as described in McCullagh (2009).

`generalized_wishart` and
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
share the same \\\Sigma\\ parameterization but differ in what `S`
represents: this function takes a **squared-distance** matrix;
`wishart_covariance` takes a **covariance** matrix. When the covariance
and its implied distance representation carry the same information under
the same centering, scaling, sites, and `nu`, the two functions return
the same objective (up to numerical tolerance) because both use the same
site-contrast likelihood and omit response-only normalizing terms. This
equality does not make AIC values from arbitrary distance and covariance
fits interchangeable.

\\F\_{ST}\\ ratio estimators do not have Wishart degrees of freedom. A
passed geometry check is necessary but does not establish this sampling
likelihood. Prefer covariance-derived squared distances and document the
construction; use
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) for
ratio-estimator responses.

The function checks `S` with
[`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md)
before evaluating the likelihood. Substantive violations stop the fit.
Tiny structural deviations within the numerical tolerance are
symmetrized and set to zero where appropriate. No Euclidean correction
is applied.

## Numerical nugget limits

Fitting uses the same small-nugget safeguards as
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
retaining supplied bounds and fixed parameters. A verified numerical
zero residual-variance limit is recorded separately from a zero graph
contribution;
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
omits Wald uncertainty for its log-scale nugget coordinate.

## Effective Wishart information

For Wishart likelihoods, `nu` scales the log likelihood. Point estimates
do not depend on it, but standard errors, likelihood-ratio tests, and
information criteria do. It is a dispersion parameter that you set,
rather than the number of markers. Linkage and the population history
shared by loci can make effective information much smaller than the SNP
count. In unpublished forward-time simulations run during development,
standard errors at the SNP count were three to five times smaller than
the spread among replicate histories. Those model-based errors
understated uncertainty.

Use
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
to compare conductance formulas, and
[`terradish_rescale_nu`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)
to report inference across plausible values. For the same successful
folds and nuisance settings, changing `nu` rescales Wishart CV scores
without changing their ranking. Neither procedure identifies the true
effective information or repairs an unsuitable response matrix.

## References

McCullagh P. 2009. Marginal likelihood for distance matrices. Statistica
Sinica 19:631-649.

## See also

[`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
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

# Generate an admissible squared-distance response from a positive
# definite covariance, then evaluate the likelihood at known values.
set.seed(1)
Sigma <- 0.8 * laplacian_inv + 0.2 * diag(nrow(laplacian_inv))
observed_covariance <- rWishart(1, df = 40, Sigma = Sigma)[, , 1] / 40
observed_squared_distance <- dist_from_cov(observed_covariance)
check_distance_response(observed_squared_distance)$admissible
#> [1] TRUE
value <- generalized_wishart(laplacian_inv, observed_squared_distance,
                             nu = 40,
                             phi = c(tau = 0.8, sigma = log(0.2)))
value$objective  # negative log-likelihood; sigma is the LOG nugget
#> [1] -23.55683
```
