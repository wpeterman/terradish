# Legacy radish grid wrapper

Deprecated compatibility wrapper retained for older code that still
calls `radish_grid()`.

## Usage

``` r
radish_grid(...)
```

## Arguments

- ...:

  Arguments passed through the deprecated `radish_grid()` compatibility
  wrapper to
  [`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md).

## Value

A list of class `terradish_grid` (also inheriting `radish_grid`)
containing the parameter grid `theta`, the maximized log-likelihood
`loglik` and nuisance parameters `phi` at each grid row, and optionally
the focal-site covariance array `covariance`. Larger log-likelihoods
indicate better fit within the supplied grid. See
[`terradish_grid`](https://wpeterman.github.io/terradish/reference/terradish_grid.md).
The `approximation` element records any coarse-raster screening and
exact refinement used in the evaluation.

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
