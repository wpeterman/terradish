# Recommend optimizer and solver settings for a data set

Profiles a `terradish` graph and, optionally, runs short probe
benchmarks to recommend optimizer, line-search, solver, and
approximation settings for a fitted model. Use this function before a
large or slow fit to identify the fastest computational configuration
for your data.

## Usage

``` r
terradish_assess_settings(
  formula,
  data,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  nu = NULL,
  theta = NULL,
  nonnegative = TRUE,
  optimizer_probe = TRUE,
  solver_probe = TRUE,
  coarse_probe = FALSE,
  optimizer_candidates = c("newton_hager_zhang", "bfgs_hager_zhang", "bfgs_armijo"),
  probe_maxit = 2L,
  factorization = c("simplicial_ldl", "simplicial_ll", "supernodal_ll"),
  solve_backend = c("matrix", "cholmod_cpp", "cholmod_cpp_cached"),
  coarse_factor = 2L,
  coarse_refine_maxit = 1L,
  loglik_tolerance = 0.001,
  max_seconds = 120,
  cores = 1L,
  verbose = FALSE
)
```

## Arguments

- formula:

  Model formula passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

- data:

  A `terradish_graph` object returned by
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md).

- conductance_model:

  Conductance-model factory (e.g.
  [`loglinear_conductance`](https://wpeterman.github.io/terradish/reference/loglinear_conductance.md)).

- measurement_model:

  Measurement model used during optimizer probes. Should match the model
  you plan to use for the final fit.

- nu:

  Effective Wishart degrees of freedom passed to measurement models that
  require it
  ([`generalized_wishart`](https://wpeterman.github.io/terradish/reference/generalized_wishart.md),
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md)).
  It is a dispersion parameter that you set, rather than the number of
  markers: linkage and the population history shared by loci can make
  effective information much smaller than the SNP count. Use the same
  value you plan to use for the final fit; it does not change the speed
  comparison. See the effective-information section of
  [`wishart_covariance`](https://wpeterman.github.io/terradish/reference/wishart_covariance.md).
  Ignored by
  [`mlpe`](https://wpeterman.github.io/terradish/reference/mlpe.md) and
  [`leastsquares`](https://wpeterman.github.io/terradish/reference/leastsquares.md).

- theta:

  Optional starting conductance parameter values. If `NULL`, the default
  all-zero start is used.

- nonnegative:

  Should regression measurement models constrain the IBR slope to be
  nonnegative during probes?

- optimizer_probe:

  Logical. Run short optimization probes to compare Newton-Raphson vs.
  BFGS and line-search algorithms?

- solver_probe:

  Logical. Benchmark direct Cholesky solver settings (factorization
  style and solve backend) via
  [`terradish_solver_benchmark`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)?

- coarse_probe:

  Logical. Include a coarse-raster warm-start probe? Only applicable
  when `data` retains its raster stack (`saveStack = TRUE`).

- optimizer_candidates:

  Character vector of optimizer configurations to probe. Supported
  entries: `"newton_hager_zhang"`, `"bfgs_hager_zhang"`,
  `"bfgs_armijo"`.

- probe_maxit:

  Maximum optimizer steps per probe. Keep small (default `2`); the
  assessment compares early iteration cost, not convergence.

- factorization:

  Character vector of Cholesky factorization styles to benchmark when
  `solver_probe = TRUE`. Options: `"simplicial_ldl"`, `"simplicial_ll"`,
  `"supernodal_ll"`.

- solve_backend:

  Character vector of solve backends to benchmark when
  `solver_probe = TRUE`. Options: `"matrix"`, `"cholmod_cpp"`,
  `"cholmod_cpp_cached"`.

- coarse_factor:

  Raster aggregation factor for the coarse-raster probe (passed to the
  `factor` entry of `approximation_control`).

- coarse_refine_maxit:

  Maximum full-resolution refinement steps for the coarse-raster probe.

- loglik_tolerance:

  Maximum log-likelihood difference below the best probe that is still
  considered acceptable when choosing the fastest option.

- max_seconds:

  Soft wall-time budget in seconds. Each probe that starts before the
  budget is exceeded is allowed to finish; later probes are skipped.

- cores:

  Number of parallel worker processes passed to probes.

- verbose:

  Logical. Print progress as probes run?

## Value

An object of class `"terradish_setting_assessment"` with components:

- `profile`:

  A list of graph statistics (number of vertices, edges, focal sites,
  parameters, and, for
  [`smooth_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/smooth_loglinear_conductance.md),
  spline-expansion counts.

- `benchmarks`:

  A list of probe results: `optimizer` (a data frame with one row per
  optimizer candidate, or `NULL` if `optimizer_probe = FALSE`), `solver`
  (a data frame of solver probe results, or `NULL` if
  `solver_probe = FALSE`), and `coarse` (the coarse-raster probe, or
  `NULL`).

- `reasons`:

  The reason recorded for each recommended setting.

- `call`:

  The matched call.

- `defaults`:

  The terradish defaults for the profiled model, after resolving
  `optimizer = "auto"` on the expanded conductance parameter count.

- `comparison`:

  A compact comparison between the recommended settings and the
  terradish defaults, including a one-line summary.

- `recommended`:

  A named list of recommended settings: `optimizer`, `control`,
  `solver`, `solver_control`, `approximation`, `approximation_control`.

- `notes`:

  Character vector of advisory messages.

A `print` method summarizes the profile and recommendations.

## Details

`terradish_assess_settings` is an advisory tool, not a model-selection
procedure. It does not change package defaults and its output does not
affect parameter estimation. Use it to choose *computational* settings;
then fit and interpret your model as usual. It is a speed probe. A fast
candidate is not evidence of optimizer convergence, identification,
correct standard errors, or biological validity.

Short probes work well for comparing early-iteration cost of different
optimizers. They are most informative for large graphs (\> 100,000
vertices) where a single function evaluation is slow. For small graphs
the default Newton-Raphson + Hager-Zhang settings typically converge in
seconds, and assessment adds more overhead than it saves.

The returned `$recommended` list is intentionally a plain named list of
arguments (`optimizer`, `control`, `solver`, `solver_control`,
`approximation`, `approximation_control`). You can inspect, edit, or
pass it directly:

    fit <- terradish(formula, data,
                     optimizer = assessment$recommended$optimizer,
                     control   = assessment$recommended$control,
                     solver    = assessment$recommended$solver)

## See also

[`terradish_solver_benchmark`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md),
[`slim_terradish`](https://wpeterman.github.io/terradish/reference/slim_terradish.md)

Other solver settings:
[`print.terradish_setting_assessment()`](https://wpeterman.github.io/terradish/reference/print.terradish_setting_assessment.md),
[`slim_terradish()`](https://wpeterman.github.io/terradish/reference/slim_terradish.md),
[`terradish_solver_benchmark()`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)

## Examples

``` r
# \donttest{
# Coarsened for the example. Assessment is most useful on the large graphs
# where solver choice actually matters, so use the full resolution there.
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

assessment <- terradish_assess_settings(
  Fst ~ forestcover + altitude,
  data = surface,
  measurement_model = mlpe,
  probe_maxit = 2
)
assessment
#> terradish setting assessment
#> Graph: 2531 vertices, 9815 edges, 31 focal sites
#> Conductance parameters: 2 
#> 
#> Recommended settings:
#>   optimizer: newton (fastest acceptable probe: newton_hager_zhang)
#>   line search: HagerZhangControl
#>   solver: direct (fastest direct solver probe: simplicial_ll + cholmod_cpp_cached)
#>   solver_control:
#> $factorization
#> [1] "simplicial_ll"
#> 
#> $solve_backend
#> [1] "cholmod_cpp_cached"
#> 
#>   approximation: none (exact full-resolution fit)
#> 
#> Defaults comparison: Assessment differs from terradish defaults: solver_control tuned from terradish defaults.
#> 
#> Optimizer probes:
#>               label elapsed loglik gradient_evaluations hessian_evaluations
#>  newton_hager_zhang    0.22   1495                   10                   4
#>    bfgs_hager_zhang    0.13   1494                    8                   1
#>         bfgs_armijo    0.09   1492                    5                   1
#>  line_search_trials status
#>                   6     OK
#>                   5     OK
#>                   2     OK
#> 
#> Fastest direct solver probes:
#>   factorization      solve_backend setup_time solve_time total_time
#>   simplicial_ll cholmod_cpp_cached          0  0.0009573       0.00
#>   supernodal_ll cholmod_cpp_cached          0  0.0035267       0.00
#>  simplicial_ldl        cholmod_cpp          0  0.0011158       0.01
#> 

fit <- terradish(
  Fst ~ forestcover + altitude,
  data = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer = assessment$recommended$optimizer,
  control = assessment$recommended$control,
  solver = assessment$recommended$solver,
  solver_control = assessment$recommended$solver_control
)

smooth_assessment <- terradish_assess_settings(
  Fst ~ forestcover + s(altitude, df = 3),
  data = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  probe_maxit = 1
)
smooth_assessment$comparison$summary
#> [1] "Assessment differs from terradish defaults: optimizer `bfgs` -> `newton`; solver_control tuned from terradish defaults."
# }
```
