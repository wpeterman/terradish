# Legacy terradish scale-optimization alias

Deprecated compatibility wrapper retained for older code that still
calls `terra_radish_scale_optim()`. It signals a deprecation warning and
forwards its arguments unchanged.

## Usage

``` r
terra_radish_scale_optim(...)
```

## Arguments

- ...:

  Arguments passed through the deprecated `terra_radish_scale_optim()`
  compatibility wrapper to
  [`terradish_scale_optim`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md).

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

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
