# Simulate genetic covariance matrices from a known conductance surface

Simulates one or more covariance response matrices from a known
conductance-resistance relationship, using the same Wishart covariance
model fitted by
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).

## Usage

``` r
simulate_covariance_response(
  theta,
  formula,
  data,
  conductance_model = loglinear_conductance,
  tau = 1,
  sigma = 0,
  nu,
  nsim = 1,
  seed = NULL,
  cores = 1L
)
```

## Arguments

- theta:

  Conductance parameters. May be supplied as a numeric vector or a
  one-row matrix. If unnamed, values are matched in the order implied by
  `formula`. Gaussian `sigma.<layer>` parameters use the map units
  reported by [`coef()`](https://rdrr.io/r/stats/coef.html); conversion
  to optimizer units is internal.

- formula:

  Model formula specifying the conductance covariates. The left hand
  side, if supplied, is ignored.

- data:

  A
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
  object.

- conductance_model:

  Conductance-model factory, such as
  [`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md).

- tau:

  Nonnegative scaling applied to the conductance-implied covariance
  matrix, on the same scale as the fitted Wishart `tau`.

- sigma:

  Nonnegative nugget variance added to the diagonal, on its natural
  (variance) scale. This differs from the fitted `sigma` reported for
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
  and
  [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
  which is the logarithm of the nugget variance. To simulate from a
  fitted model, pass [`exp()`](https://rdrr.io/r/base/Log.html) of the
  fitted `sigma`.

- nu:

  Effective Wishart degrees of freedom: controls the amount of sampling
  noise in the simulated covariance: larger `nu` gives tighter draws
  around `Sigma`. This is a model assumption, not a marker count or an
  estimate of information from biological data. Examine a range of
  plausible values. The simulation requires `nu` to be at least the
  covariance-matrix dimension. See
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
  for how `nu` propagates into estimation (it scales standard errors and
  model-selection statistics but not point estimates).

- nsim:

  Number of covariance matrices to simulate.

- seed:

  Optional random seed.

- cores:

  Number of cores passed to
  [`terradish_distance`](https://wpeterman.github.io/terradish/reference/terradish_distance.md).

## Value

A list containing:

- covariance:

  A covariance matrix if `nsim = 1`, otherwise a three-dimensional array
  of simulated covariance matrices.

- Sigma:

  The covariance matrix used as the Wishart scale matrix after applying
  `tau` and `sigma`.

- E:

  The conductance-implied covariance matrix returned by
  [`terradish_distance`](https://wpeterman.github.io/terradish/reference/terradish_distance.md)`(covariance = TRUE)`.

- theta:

  The validated conductance parameter matrix used for simulation.

- tau:

  The supplied `tau`.

- sigma:

  The supplied `sigma` (nugget variance, natural scale).

- nu:

  The supplied `nu`.

## Details

The helper first computes the conductance-implied covariance matrix
`E(theta)` from the supplied resistance surface. It then forms
`Sigma = tau * E(theta) + sigma * I` and simulates
`S ~ Wishart(nu, Sigma) / nu`. This makes the returned covariance
matrices directly compatible with
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).
These results are internal to the assumed model. Simulating and fitting
at the same incorrect `nu` cannot reveal that error. Use
[`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md)
with different `nu` and `nu_fit` to examine its effect on coverage.

## See also

[`covariance_response_power`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)

Other simulation and study design:
[`covariance_response_power()`](https://wpeterman.github.io/terradish/reference/covariance_response_power.md),
[`print.terradish_covariance_power()`](https://wpeterman.github.io/terradish/reference/print.terradish_covariance_power.md)

## Examples

``` r
r1 <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0,
                  ymax = 3, vals = 1:9)
r2 <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0,
                  ymax = 3, vals = c(9, 1, 4, 3, 8, 2, 5, 7, 6))
covariates <- c(r1, r2)
names(covariates) <- c("x1", "x2")
pts <- terra::vect(matrix(c(0.5, 0.5,
                           1.5, 1.5,
                           2.5, 2.5), ncol = 2, byrow = TRUE),
                   type = "points")
surface <- conductance_surface(covariates, pts, directions = 4)
sim <- simulate_covariance_response(
  theta = c(x1 = 0.3, x2 = -0.2),
  formula = ~ x1 + x2,
  data = surface,
  tau = 0.8,
  sigma = 0.5,
  nu = 20,
  seed = 1
)
sim$covariance
#>            [,1]       [,2]       [,3]
#> [1,]  0.4538078  0.1349725 -0.2132605
#> [2,]  0.1349725  0.7820864 -0.2164368
#> [3,] -0.2132605 -0.2164368  0.7259222
```
