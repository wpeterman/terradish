# Legacy radish algorithm wrapper

Deprecated compatibility wrapper retained for legacy code that still
calls `radish_algorithm()` directly.

## Usage

``` r
radish_algorithm(...)
```

## Arguments

- ...:

  Arguments passed through the deprecated `radish_algorithm()`
  compatibility wrapper to
  [`terradish_algorithm`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md).

## Value

A list containing at a minimum:

- covariance:

  rows/columns of the generalized inverse of the graph Laplacian for a
  subset of target vertices

Additionally, if 'objective == TRUE':

- objective:

  (if 'objective') the negative loglikelihood

- phi:

  (if 'objective') fitted values of the nuisance parameters of 'g'

- subproblem:

  (if 'objective') convergence code and iteration count for the
  nuisance-parameter profile

- boundary:

  (if 'objective') is the solution on the boundary (e.g. no genetic
  structure)?

- fitted:

  (if 'objective') matrix of expected genetic distances among target
  vertices

- gradient:

  (if 'gradient') gradient of negative loglikelihood with respect to
  theta

- hessian:

  (if 'hessian') Hessian matrix of the negative loglikelihood with
  respect to theta

- partial_X:

  (if 'partial') Jacobian of the gradient with respect to the spatial
  covariates

- partial_S:

  (if 'partial') Jacobian of the gradient with respect to the observed
  genetic distances
