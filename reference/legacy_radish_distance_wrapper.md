# Legacy radish distance wrapper

Deprecated compatibility wrapper retained for older code that still
calls `radish_distance()`.

## Usage

``` r
radish_distance(...)
```

## Arguments

- ...:

  Arguments passed through the deprecated `radish_distance()`
  compatibility wrapper to
  [`terradish_distance`](https://wpeterman.github.io/terradish/reference/terradish_distance.md).

## Value

A list of class `terradish_grid` (also inheriting `radish_grid`)
containing the parameter grid `theta` and either `distance`, an array of
resistance distances, or `covariance`, an array of graph covariances
when `covariance = TRUE`. The first two array dimensions index focal
sites; the third indexes parameter-grid rows. See
[`terradish_distance`](https://wpeterman.github.io/terradish/reference/terradish_distance.md).

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
