# Legacy terradish scale-optimization alias

Legacy alias retained for backward compatibility with older
`terra_radish_scale_optim()` calls.

## Usage

``` r
terra_radish_scale_optim(
  formula,
  covariates,
  coords,
  scales = NULL,
  lower = NULL,
  upper = NULL,
  scale_grid = NULL,
  search = c("coordinate", "grid"),
  objective = c("aic", "logLik"),
  grid_points = 5L,
  maxit = 10L,
  tol = 0.1,
  scale_fun = "terradish",
  scale_args = NULL,
  postprocess = scale_covariates,
  directions = 8L,
  theta = NULL,
  verbose = FALSE,
  ...
)
```

## Value

An object of class `"terradish_scale_optim"`, a list with:

- par:

  The selected scale for each raster layer, in the units `scale_fun`
  uses (map units for both built-in methods).

- value:

  The outer objective at `par`. Lower is better, whichever objective was
  used.

- fit:

  The fitted `terradish` model at `par`. Read its conductance
  coefficients with [`coef()`](https://rdrr.io/r/stats/coef.html) the
  same way as any other fit; they are on the log-conductance scale.

- surface:

  The conductance surface built from the rasters smoothed at `par`,
  ready to pass to
  [`conductance()`](https://wpeterman.github.io/terradish/reference/conductance.md)
  or [`plot()`](https://rspatial.github.io/terra/reference/plot.html).

- evaluations:

  A data frame with one row per candidate scale combination: the scale
  for each layer, the `objective` value, `aic`, `loglik`, and a `status`
  column that records `"OK"` or the error for candidates that failed to
  fit. Inspect this to see how sharply the objective distinguishes the
  scales; a flat column means the objective does not identify a
  smoothing scale within the searched range.

- search:

  The search strategy used, `"grid"` or `"coordinate"`.

- objective:

  The name of the objective minimized, `"aic"` or `"logLik"`.

- iterations, convergence:

  Number of outer iterations, and `0` if the coordinate search converged
  (always `0` for a grid search).
