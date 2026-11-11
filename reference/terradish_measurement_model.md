# Measurement models: the likelihood interface

A measurement model links the graph kernel implied by a conductance
surface to the observed genetic response. It is a function of class
`"terradish_measurement_model"` that evaluates a negative
log-likelihood, and its derivatives, for a matrix of observed pairwise
genetic distances or genetic covariance given a submatrix of the
generalized inverse of the graph Laplacian. Pass one to the
`measurement_model` argument of
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Details

**Choosing a model.**

- [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md):

  Pairwise genetic distances (for example \\F\_{ST}\\), with correlation
  among pairs that share a site. The default.

- [`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md):

  Pairwise genetic distances with independent errors. Fast; for
  exploration and diagnostics.

- [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md):

  An allele-frequency covariance matrix. Requires the effective degrees
  of freedom `nu`.

- [`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md):

  Squared distances derived from a covariance matrix. Requires `nu` and
  an admissible response (see
  [`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md)).

- [`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
  [`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md):

  Constructors that return a measurement model extended with pairwise
  environmental covariates (isolation by environment).

**Interface.** Every measurement model has the signature
`g(E, S, phi, nu, gradient, hessian, partial, nonnegative, validate)`.
`E` is the focal-site block of the generalized inverse of the graph
Laplacian, `S` is the observed response matrix, `phi` is the vector of
nuisance parameters, and `nu` is the effective Wishart degrees of
freedom (unused by the distance models). Called without `phi`, the
function returns a list with starting values `phi` and the bounds
`lower` and `upper`. Called with `phi`, it returns a list with the
negative log-likelihood `objective`, the expected response `fitted`, a
`boundary` flag that is `TRUE` when the resistance coefficient is zero,
and, on request, `gradient` and `hessian` with respect to `phi` and the
derivative pieces (`gradient_E`, `partial_E`, `partial_S`, `jacobian_E`,
`jacobian_S`) that
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
uses to differentiate through the graph. The legacy class
`"radish_measurement_model"` is also carried for compatibility.

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
profiles `phi` at every value of the conductance parameters, so `phi`
are nuisance parameters. They are reported by
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
on each model's own scale; the same name can mean different things in
different models (see the parameter glossary in
[`terradish-package`](https://wpeterman.github.io/terradish/reference/terradish-package.md)).

**Wishart nuisance fitting.** For the ordinary two-component
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)
and
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md)
models with their standard bounds, nuisance fitting first profiles the
overall covariance variance analytically and searches the relative
graph/nugget contribution on site contrasts. This deterministic
initialization avoids false convergence when a warm log nugget has an
almost zero derivative. The existing Newton optimizer then refines the
result and computes derivatives at the same parameters. Fixed predictive
parameters and custom bounds retain their supplied values and
constraints. The finite profile search does not guarantee a global
optimum for every possible covariance matrix. A nugget approaching zero
remains a limiting estimate and can imply weak nuisance identification.

## Worked examples

[`vignette("v02-wishart-covariance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v02-wishart-covariance.md)
for choosing and checking a covariance response, and
[`vignette("v03-mlpe-pairwise-distance", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v03-mlpe-pairwise-distance.md)
for pairwise genetic distances.

## See also

[`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`mlpe_covariates`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`wishart_covariates`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md),
[`check_distance_response`](https://wpeterman.github.io/terradish/reference/check_distance_response.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other measurement models:
[`generalized_wishart()`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
[`leastsquares()`](https://wpeterman.github.io/terradish/reference/leastsquares.md),
[`mlpe()`](https://wpeterman.github.io/terradish/reference/mlpe.md),
[`mlpe_covariates()`](https://wpeterman.github.io/terradish/reference/mlpe_covariates.md),
[`pairwise_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_covariates.md),
[`pairwise_endpoint_covariates()`](https://wpeterman.github.io/terradish/reference/pairwise_endpoint_covariates.md),
[`wishart_covariance()`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md),
[`wishart_covariates()`](https://wpeterman.github.io/terradish/reference/wishart_covariates.md)
