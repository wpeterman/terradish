# Print and read a settings assessment

Compact display for the object returned by
[`terradish_assess_settings`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md):
the graph profile, the recommended computational settings and why each
was chosen, the probe timings behind them, and any warnings the
assessment raised.

## Usage

``` r
# S3 method for class 'terradish_setting_assessment'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `terradish_setting_assessment`, from
  [`terradish_assess_settings`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md).

- digits:

  Number of significant digits used when formatting the probe timing
  tables.

- ...:

  Ignored, present for compatibility with the `print` generic.

## Value

`x`, invisibly. Called for the printed output. The machine readable
pieces are `x$profile`, `x$recommended`, `x$reasons`, `x$benchmarks`,
`x$comparison`, and `x$notes`.

## Details

**How to read the output.** It comes in four parts.

*Graph* reports the size of the problem: vertices (active raster cells),
edges, focal sites, and conductance parameters. This is the number that
decides whether any of the rest matters. Below roughly 50,000 vertices
the defaults are almost always fine; above a few hundred thousand,
solver choice dominates the runtime.

*Recommended settings* lists the optimizer, line search, solver, and
approximation the probes favored, each followed by a parenthesized
reason. Pass them straight through with
`fit <- terradish(..., optimizer = rec$optimizer, control = rec$control, solver = rec$solver)`,
where `rec <- assessment$recommended`. These choices change only the
computational path, never the estimates.

*Defaults comparison* is the one-line verdict: whether the assessment
agreed with the package defaults. If it says it matched, there is
nothing to do and you can drop the assessment from your script.

*Probes* are the measurements behind the recommendation. The optimizer
table gives elapsed time and the log-likelihood each candidate reached
in the probe budget, so a candidate that is fast but reaches a worse
log-likelihood is not actually better. The solver table shows the three
fastest direct factorizations split into setup and solve time; setup
dominates when a fit re-factorizes often, solve time when the focal set
is large. *Notes* collects anything that needs attention, such as a
solver that failed to converge during probing.

## See also

[`terradish_assess_settings`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md)
for the assessment itself and what each component holds,
[`terradish_solver_benchmark`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)
for benchmarking solvers alone, and
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
for the fit the recommendations feed. See
[`vignette("v11-large-landscapes", package = "terradish")`](https://wpeterman.github.io/terradish/articles/v11-large-landscapes.md)
for when this is worth running.

Other solver settings:
[`slim_terradish()`](https://wpeterman.github.io/terradish/reference/slim_terradish.md),
[`terradish_assess_settings()`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md),
[`terradish_solver_benchmark()`](https://wpeterman.github.io/terradish/reference/terradish_solver_benchmark.md)

## Examples

``` r
# \donttest{
data(melip)
covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")
# coarsened for the example; assessment earns its keep on large graphs
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
coords <- terra::unwrap(melip.coords)
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords[keep])

assessment <- terradish_assess_settings(
  Fst ~ forestcover + altitude,
  data = surface, measurement_model = mlpe, probe_maxit = 2)

assessment                      # the four-part report described above
#> terradish setting assessment
#> Graph: 2531 vertices, 9815 edges, 31 focal sites
#> Conductance parameters: 2 
#> 
#> Recommended settings:
#>   optimizer: newton (fastest acceptable probe: newton_hager_zhang)
#>   line search: HagerZhangControl
#>   solver: direct (fastest direct solver probe: simplicial_ldl + cholmod_cpp)
#>   solver_control:
#> $factorization
#> [1] "simplicial_ldl"
#> 
#> $solve_backend
#> [1] "cholmod_cpp"
#> 
#>   approximation: none (exact full-resolution fit)
#> 
#> Defaults comparison: Assessment differs from terradish defaults: solver_control tuned from terradish defaults.
#> 
#> Optimizer probes:
#>               label elapsed loglik gradient_evaluations hessian_evaluations
#>  newton_hager_zhang    0.22   1495                   10                   4
#>    bfgs_hager_zhang    0.19   1494                    8                   1
#>         bfgs_armijo    0.10   1492                    5                   1
#>  line_search_trials status
#>                   6     OK
#>                   5     OK
#>                   2     OK
#> 
#> Fastest direct solver probes:
#>   factorization solve_backend setup_time solve_time total_time
#>  simplicial_ldl   cholmod_cpp          0  0.0008286          0
#>   simplicial_ll        matrix          0  0.0000000          0
#>   supernodal_ll        matrix          0  0.0000000          0
#> 
assessment$comparison$summary   # one-line verdict against the defaults
#> [1] "Assessment differs from terradish defaults: solver_control tuned from terradish defaults."
assessment$recommended$solver   # the settings, ready to reuse
#> [1] "direct"
# }
```
