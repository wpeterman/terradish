# Low-level likelihood and derivatives at one parameter value

Calculates likelihood, gradient, Hessian, and partial derivatives of a
parameterized conductance surface, given a function mapping spatial data
to conductance and a function mapping resistance distance (covariance)
to the genetic response. This is the evaluation that
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
calls at every optimizer step; most users never call it directly. The
approach of fitting one parameterized surface is described by Peterman
and Pope (2021).

## Usage

``` r
terradish_algorithm(
  f,
  g,
  s,
  S,
  theta,
  nu = NULL,
  phi = NULL,
  objective = TRUE,
  gradient = TRUE,
  hessian = TRUE,
  partial = TRUE,
  nonnegative = TRUE,
  validate = FALSE,
  cores = 1L,
  curvature = c("exact", "gauss_newton"),
  solver = c("direct", "auto", "amg"),
  solver_control = NULL,
  solver_warm_start = NULL,
  solver_reuse_state = NULL,
  measurement_control = NULL,
  worker_pool = NULL
)
```

## Arguments

- f:

  A function of class `"terradish_conductance_model"`, as returned by a
  conductance-model factory such as
  [`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md)
  (see
  [`terradish_conductance_model_factory`](https://wpeterman.github.io/terradish/reference/terradish_conductance_model_factory.md)).

- g:

  A function of class `"terradish_measurement_model"` (see
  [`terradish_measurement_model`](https://wpeterman.github.io/terradish/reference/terradish_measurement_model.md)).

- s:

  An object of class `"terradish_graph"`

- S:

  A matrix of observed genetic distances or covariance, as required by
  `g`

- theta:

  Parameters for conductance surface (e.g. inputs to 'f')

- nu:

  Effective Wishart degrees of freedom, passed through to measurement
  models that require it.

- phi:

  Optional warm-start values for the nuisance-parameter subproblem.

- objective:

  Compute negative loglikelihood?

- gradient:

  Compute gradient of negative loglikelihood wrt theta?

- hessian:

  Compute Hessian matrix of negative loglikelihood wrt theta?

- partial:

  Compute partial derivatives of negative loglikelihood wrt theta and
  spatial covariates/observed genetic distances

- nonnegative:

  Force regression-like 'measurement_model' to have nonnegative slope?

- validate:

  Numerical validation via 'numDeriv' (very slow, use for debugging
  small examples)

- cores:

  Number of worker processes to use for per-parameter derivative
  calculations. `1` evaluates serially.

- curvature:

  Which curvature to return in `hessian`. `"exact"` (default) returns
  the exact Hessian of the negative log-likelihood. `"gauss_newton"`
  returns the Gauss-Newton / Fisher-information approximation, obtained
  by dropping the two residual-weighted second-derivative terms (the
  second derivative of the resistance covariance \\E\\ and of
  conductance with respect to `theta`). It is an information-based
  approximation and is positive semidefinite when the retained
  measurement-model curvature is positive semidefinite. It equals the
  exact Hessian only when the omitted terms are zero or negligible;
  convergence alone does not guarantee equality. A large gap can reflect
  residual structure, model mismatch, or numerical instability and
  should be investigated rather than assigned a unique cause. It also
  requires only first derivatives of the conductance model, which is
  useful for Gaussian scale-aware and spline conductance models whose
  second derivatives are expensive or unstable. The number of linear
  solves is the same as for the exact Hessian.

- solver:

  Linear-system solver. `"direct"` uses sparse Cholesky, `"amg"` uses
  algebraic multigrid, and `"auto"` chooses between them.

- solver_control:

  Optional named list of solver settings. For `solver = "direct"`,
  supported entries include `factorization` (`"auto"`,
  `"simplicial_ldl"`, `"simplicial_ll"`, or `"supernodal_ll"`),
  `solve_backend` (`"matrix"` or the experimental `"cholmod_cpp"` and
  `"cholmod_cpp_cached"` backends), `supernodal_min_vertices`,
  `supernodal_max_rhs`, and `perm`. For `solver = "auto"`, supported
  selection entries include `auto_direct_max_vertices`,
  `auto_amg_min_vertices`. The legacy `auto_direct_max_rhs` is ignored.
  For `solver = "amg"`, supported entries include `tol`, `maxit`,
  `coarse_enough`, `npre`, `npost`, `sa_relax`, `aggr_eps_strong`,
  `estimate_spectral_radius`, `power_iters`, and `reuse_preconditioner`.
  `reuse_preconditioner_max_age` can be set to a finite nonnegative
  value to periodically rebuild the AMG hierarchy instead of reusing it
  indefinitely. Direct solvers also accept `reference = "auto"`
  (default), which grounds a maximum-conductance vertex when the
  conductance ratio exceeds one million, or `"last"` and
  `"max_conductance"`. This changes the numerical reference potential
  while preserving the graph and model. It improves some stiff solves
  but cannot guarantee accuracy at arbitrary conductance contrasts.
  Direct solves check the maximum absolute residual in each right-hand
  side, divided by its maximum absolute input or one, whichever is
  larger. `residual_tolerance = 1e-8` rejects unresolved systems by
  default; `Inf` disables this check for diagnosis. A small residual is
  necessary but does not guarantee small forward error on an
  ill-conditioned graph. Before rejection, `refinement_steps = 3`
  permits up to three correction solves using the same factorization and
  the current residual. Direct supernodal factorizations can benefit
  from a threaded BLAS, but the relevant thread counts are controlled by
  the external R/BLAS build rather than by `terradish_algorithm()`. If
  the direct solver stops with a CHOLMOD or non-positive-definite
  reduced Laplacian message, the current conductance values made the
  graph numerically singular. Check for missing or infinite raster
  covariates, use scaled covariates and moderate starting values near
  zero, and try `solver = "auto"` or `solver = "amg"` while diagnosing
  the fit.

- solver_warm_start:

  Optional initial guess for the reduced-system solve. This is primarily
  useful for iterative solvers when evaluating nearby parameter values.

- solver_reuse_state:

  Optional reusable solver state returned by a prior
  `terradish_algorithm()` call. This is used to reuse AMG hierarchy
  information or compatible direct CHOLMOD factorization state across
  nearby evaluations.

- measurement_control:

  Optional
  [`NewtonRaphsonControl`](https://wpeterman.github.io/terradish/reference/NewtonRaphsonControl.md)
  object used to profile nuisance parameters. The default retains the
  package's high-accuracy profiling tolerances.

- worker_pool:

  Optional reusable worker pool created internally by
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).
  Most users should leave this as `NULL`.

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

## References

Peterman WE, Pope NS (2021). The use and misuse of regression models in
landscape genetic analyses. Molecular Ecology, 30, 37–47.
[doi:10.1111/mec.15716](https://doi.org/10.1111/mec.15716) .

## See also

[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md)

Other model fitting:
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md),
[`terradish()`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`terradish_distance()`](https://wpeterman.github.io/terradish/reference/terradish_distance.md),
[`terradish_grid()`](https://wpeterman.github.io/terradish/reference/terradish_grid.md),
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
conductance_model <- loglinear_conductance(~ altitude + forestcover, surface$x)

# One evaluation at theta = (-0.3, 0.3): negative log-likelihood, its
# gradient and Hessian with respect to theta, and the profiled nuisance
# parameters.
value <- terradish_algorithm(conductance_model, leastsquares,
                             surface, Fst, theta = c(-0.3, 0.3))
value$objective
#> [1] -1375.273
value$gradient
#> [1] -94.60471  51.19362
value$phi
#>              [,1]
#> alpha -0.01942733
#> beta   0.31180555
#> tau    6.91515095
# }
```
