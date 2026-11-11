# Control settings for the Newton and BFGS optimizers

Tuning parameters shared by Newton-Raphson and BFGS quasi-Newton
optimizers inside
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Usage

``` r
NewtonRaphsonControl(
  maxit = 100,
  ctol = sqrt(.Machine$double.eps),
  ftol = sqrt(.Machine$double.eps),
  etol = 10 * .Machine$double.eps,
  verbose = FALSE,
  eps = 1e-08,
  del = 1,
  ls.control = HagerZhangControl()
)
```

## Arguments

- maxit:

  Maximum number of Newton or quasi-Newton steps. Increase if the
  optimizer reports a convergence warning on large or difficult
  problems.

- ctol:

  Gradient convergence tolerance. Optimization stops when the largest
  absolute projected gradient is below `ctol`.

- ftol:

  Objective convergence tolerance. Optimization also stops when the
  absolute change in the objective is below `ftol` and the largest
  absolute projected gradient is below `sqrt(ctol)`.

- etol:

  Eigenvalue threshold used to detect near-singular Hessians. Absolute
  eigenvalues smaller than `etol` times the largest absolute Hessian
  entry are replaced by one when forming the Newton step. Increasing
  this value regularizes the step at the cost of slower convergence near
  flat regions.

- verbose:

  Logical. If `TRUE`, report the iteration count, objective value, and
  gradient norm at each step (emitted with
  [`message()`](https://rdrr.io/r/base/message.html)). The `verbose`
  argument of
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
  overrides this field when it is set explicitly.

- eps:

  Box-constraint tolerance. A parameter is considered to be "on the
  boundary" when it is within `eps * abs(par)` of the bound.
  Active-constraint gradient components are zeroed before the step is
  taken.

- del:

  Typical initial step length for quasi-Newton (BFGS) methods. Used to
  scale the initial approximate Hessian before any updates are
  available.

- ls.control:

  Control object for the line search, created by
  [`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md)
  (default) or
  [`ArmijoControl`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md).
  `ArmijoControl` is only valid when `optimizer = "bfgs"` is passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Value

A named list holding the settings passed to the optimizer (a plain list
with no additional class). The elements are the arguments above,
validated and stored unchanged: `maxit`, `ctol`, `ftol`, `etol`,
`verbose`, `eps`, `del`, and `ls.control` (itself a `HagerZhangControl`
or `ArmijoControl` list). Pass the object to the `control` argument of
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md);
inspect it with [`str()`](https://rdrr.io/r/utils/str.html). The list
carries no fitted quantities, so nothing in it needs back-transforming.

## Details

The projected gradient sets outward-pointing components at active bounds
to zero. A sufficiently small projected gradient gives convergence code
0. A small objective change also gives code 0 if the projected gradient
is below `sqrt(ctol)`. Otherwise the optimizer warns that it has stalled
and returns code 2. Reaching the iteration limit gives code 1.

For most resistance-surface models the defaults converge reliably.
Consider changing them when:

- The optimizer hits `maxit` before converging: increase `maxit`.

- You want a quick exploratory fit: reduce `maxit` (for example 5 to
  10).

- The optimizer reports a stall (code 2) or repeated line-search
  failures: scale the covariates, start closer to zero, or adjust the
  line search through `ls.control`.

## See also

[`HagerZhangControl`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md),
[`ArmijoControl`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md),
[`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)

Other optimizer controls:
[`ArmijoControl()`](https://wpeterman.github.io/terradish/reference/ArmijoControl.md),
[`HagerZhangControl()`](https://wpeterman.github.io/terradish/reference/HagerZhangControl.md)

## Examples

``` r
ctrl <- NewtonRaphsonControl(maxit = 25, verbose = FALSE)
str(ctrl)
#> List of 8
#>  $ maxit     : num 25
#>  $ ctol      : num 1.49e-08
#>  $ etol      : num 2.22e-15
#>  $ ftol      : num 1.49e-08
#>  $ verbose   : logi FALSE
#>  $ eps       : num 1e-08
#>  $ del       : num 1
#>  $ ls.control:List of 10
#>   ..$ delta        : num 0.1
#>   ..$ sigma        : num 0.9
#>   ..$ alphamax     : num Inf
#>   ..$ rho          : num 5
#>   ..$ epsilon      : num 1e-06
#>   ..$ gamma        : num 0.66
#>   ..$ linesearchmax: num 50
#>   ..$ psi3         : num 0.1
#>   ..$ c            : num 1
#>   ..$ verbose      : logi FALSE

# Relaxed tolerances for a quick exploratory run:
ctrl_quick <- NewtonRaphsonControl(maxit = 10, ctol = 1e-3, ftol = 1e-3)

# \donttest{
data(melip)
covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
                                                na.rm = TRUE))
# Keep one site per raster cell (some melip sites share a cell).
coords <- terra::unwrap(melip.coords)
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
Fst <- melip.Fst[keep, keep]
surface <- conductance_surface(covariates, coords[keep])

# A tighter, longer-running control than the default
fit <- terradish(Fst ~ forestcover + altitude, data = surface,
                 conductance_model = loglinear_conductance,
                 measurement_model = mlpe,
                 control = NewtonRaphsonControl(maxit = 200, ctol = 1e-8))
coef(fit)
#> forestcover    altitude 
#> -0.09612631 -0.79206795 

# The same fit, reporting each optimizer step. `verbose` is the simple
# switch; it overrides whatever `control` says.
fit_loud <- terradish(Fst ~ forestcover + altitude, data = surface,
                      conductance_model = loglinear_conductance,
                      measurement_model = mlpe,
                      control = ctrl_quick, verbose = TRUE)
#> Projected Newton-Raphson with Hager-Zhang line search
#> [1] f(x) = 1.443e+03  |f(x)-fold(x)| = 0.000e+00  max|f'(x)| = 1.291e+02  |f''(x)| = 3.940e+03
#> [2] f(x) = 1.493e+03  |f(x)-fold(x)| = 4.998e+01  max|f'(x)| = 1.530e+01  |f''(x)| = -5.248e+03
#> [3] f(x) = 1.495e+03  |f(x)-fold(x)| = 1.722e+00  max|f'(x)| = 2.609e+00  |f''(x)| = -3.428e+03
#> [4] f(x) = 1.495e+03  |f(x)-fold(x)| = 8.574e-02  max|f'(x)| = 3.715e-01  |f''(x)| = -3.215e+03
#> [5] f(x) = 1.495e+03  |f(x)-fold(x)| = 9.996e-04  max|f'(x)| = 2.617e-03  |f''(x)| = -3.157e+03
#> Solution on interior with `max(abs(gradient))` == 0.00261731205736565 and `diff(f)` == 0.000999619663843987
# }
```
