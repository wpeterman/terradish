# Legacy radish plot method

Deprecated S3 method retained so legacy `"radish"` fit objects still
dispatch to
[`plot.terradish`](https://wpeterman.github.io/terradish/reference/plot.terradish.md).

## Usage

``` r
# S3 method for class 'radish'
plot(x, ...)
```

## Arguments

- x:

  A fitted legacy `radish` object.

- ...:

  Additional graphical parameters.

## Value

A plot object that is returned visibly so base R auto-printing works as
expected:

- `"fit"`: a `ggplot` object with observed-vs-fitted pairwise distances
  or genetic covariance, according to the response model.

- `"surface"`, `"marginal"`, `"marginal_response"`, and `"sigma"`:
  either a single `ggplot` object (one panel) or a named list of
  `ggplot` objects (multiple panels), with one plot per panel/covariate.

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
