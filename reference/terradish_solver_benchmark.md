# Benchmark sparse Cholesky solver settings on a graph

Times the linear-system setup and solve used internally by
[`terradish_algorithm`](https://wpeterman.github.io/terradish/reference/terradish_algorithm.md)
for a fixed graph and conductance vector. This is intended as a
lightweight hardware-specific diagnostic for choosing direct
factorization settings on large rasters.

## Usage

``` r
terradish_solver_benchmark(
  data,
  conductance = NULL,
  factorization = c("simplicial_ldl", "simplicial_ll", "supernodal_ll"),
  solve_backend = "matrix",
  n_replicates = 1L,
  perm = TRUE
)
```

## Arguments

- data:

  A `terradish_graph` object.

- conductance:

  Optional numeric conductance vector with one value per graph vertex.
  Defaults to constant conductance.

- factorization:

  Direct solver factorization modes to compare. Supported values are
  `"simplicial_ldl"`, `"simplicial_ll"`, and `"supernodal_ll"`.

- solve_backend:

  Direct-solver RHS solve backend to compare. `"matrix"` uses Matrix's
  standard [`solve()`](https://rdrr.io/r/base/solve.html) method. The
  experimental `"cholmod_cpp"` backend calls Matrix's CHOLMOD C API
  directly from compiled code after the factorization has been built.
  The experimental `"cholmod_cpp_cached"` backend also keeps the CHOLMOD
  factorization/update state inside a compiled external pointer.

- n_replicates:

  Number of times to repeat each setting.

- perm:

  Should CHOLMOD use a fill-reducing permutation?

## Value

A data frame with one row per benchmark replicate and columns for graph
size, factorization mode, solve backend, setup time, solve time, and
total time.

## Details

The benchmark isolates the sparse linear algebra step. It does not run a
full likelihood optimization, line search, measurement-model subproblem,
or Gaussian smoothing update. Use it to learn whether
`"simplicial_ldl"`, `"simplicial_ll"`, or `"supernodal_ll"` is fastest
on a particular machine and graph size, then pass the preferred choice
through `solver_control = list(factorization = "...")` in
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
when using the direct solver. Use `solve_backend` to evaluate whether
either experimental C++ CHOLMOD path is faster for the right-hand-side
count used by the graph. The cached backend is most informative when
compared inside a full optimization, where nearby evaluations can reuse
the same symbolic analysis.

Very small rasters can give noisy timings because setup and R overhead
are a large fraction of elapsed time. The comparison is most informative
on graphs large enough that Cholesky factorization is a visible part of
the fit time.

## See also

[`terradish_assess_settings`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other solver settings:
[`print.terradish_setting_assessment()`](https://wpeterman.github.io/terradish/reference/print.terradish_setting_assessment.md),
[`slim_terradish()`](https://wpeterman.github.io/terradish/reference/slim_terradish.md),
[`terradish_assess_settings()`](https://wpeterman.github.io/terradish/reference/terradish_assess_settings.md)

## Examples

``` r
library(terra)
data(melip)
altitude <- terra::unwrap(melip.altitude)
coords <- terra::unwrap(melip.coords)
names(altitude) <- "altitude"
# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(altitude, terra::crds(coords)))
surface <- conductance_surface(altitude, coords[keep])
# One row per setting: setup (factorization) and solve times in seconds.
terradish_solver_benchmark(
  surface,
  factorization = c("simplicial_ldl", "simplicial_ll"),
  n_replicates = 1
)
#>   solver  factorization solve_backend replicate n_vertices n_rhs setup_time
#> 1 direct simplicial_ldl        matrix         1      22443    34       0.05
#> 2 direct  simplicial_ll        matrix         1      22443    34       0.05
#>   solve_time total_time
#> 1       0.02       0.11
#> 2       0.02       0.11
```
