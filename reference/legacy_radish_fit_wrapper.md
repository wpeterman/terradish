# Legacy radish fit wrapper

Deprecated compatibility wrapper retained for older code that still
calls `radish()`.

## Usage

``` r
radish(...)
```

## Arguments

- ...:

  Arguments passed through the deprecated `radish()` compatibility
  wrapper to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

## Value

An object of class `c("terradish", "radish")`, a list. Most users need
only the methods:
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
for estimates, standard errors, and nuisance-parameter intervals,
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) for conductance
coefficients,
[`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
and [`plot()`](https://rspatial.github.io/terra/reference/plot.html) for
maps and curves, and [`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`anova()`](https://rdrr.io/r/stats/anova.html), or
[`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md)
for comparison (see
[`terradish_methods`](https://wpeterman.github.io/terradish/reference/terradish_methods.md)).
The list elements are:

- call, formula:

  The matched call and the model formula.

- dim:

  Named vector with the number of graph vertices (raster cells), focal
  sites, and graph edges.

- convergence:

  List describing the conductance optimization: `code` (0 converged, 1
  iteration limit reached, 2 stalled or line search failed), `message`,
  `iterations`, `max_abs_projected_gradient`, `criterion`, `boundary`
  (`TRUE` when an estimate sits on a bound or no resistance structure
  was detected), and `restarted`. Check `object$convergence$code == 0`
  before interpreting estimates.

- fit:

  The final graph and likelihood evaluation. `object$fit$phi` holds the
  nuisance (measurement-model) estimates on the measurement model's
  internal scale, for example log precision;
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  reports them with intervals. `object$fit$fitted` is the expected
  response matrix among focal sites. `object$fit$subproblem$convergence`
  is the nuisance optimizer's code (0 means converged), which should
  also be 0 before estimates are interpreted.

- loglik, df, aic:

  Maximized log-likelihood, number of estimated parameters (conductance
  coefficients plus nuisance parameters), and \\AIC = -2\ell + 2K\\.

- mle:

  Conductance estimates (`theta`, on the log-conductance scale per unit
  of each covariate as supplied to the graph), and the gradient and
  Hessian of the log-likelihood at the estimate. Elements ending in
  `_internal` are on the optimizer's scale.

- submodels:

  The conductance-model and measurement-model functions used for the
  fit;
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  and [`plot()`](https://rspatial.github.io/terra/reference/plot.html)
  reuse them.

- comparison:

  Measurement family, `nu`, graph fingerprint, and response columns that
  [`aic_table`](https://wpeterman.github.io/terradish/reference/aic_table.md)
  and [`anova()`](https://rdrr.io/r/stats/anova.html) use to check that
  fits are comparable.

- cost, diagnostics, approximation:

  Newton steps and function calls, solver timing and evaluation
  diagnostics, and a record of any landmark or coarse-raster warm start.

- gaussian_scale_info, spline_monotonicity:

  Metadata for Gaussian smoothing and spline conductance terms, or
  `NULL`.

- leverage:

  Leverage arrays for the response (`S`) and spatial covariates (`X`)
  when `leverage = TRUE`, otherwise `NULL` entries.

- storage:

  Present only when `slim = TRUE`; records the model closures and
  leverage arrays that were removed.

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
