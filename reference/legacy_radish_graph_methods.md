# Legacy radish graph conductance method

Deprecated S3 method retained so legacy `"radish_graph"` objects still
dispatch to
[`conductance`](https://wpeterman.github.io/terradish/reference/conductance.md).

## Usage

``` r
# S3 method for class 'radish_graph'
conductance(
  x,
  fit,
  quantile = 0.95,
  support = c("none", "focal"),
  support_probs = c(0, 1),
  clamp_covariates = NULL,
  ...
)
```

## Arguments

- x:

  A legacy `radish_graph` object.

- fit:

  A fitted model returned by
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).

- quantile:

  Confidence level for the pointwise conductance interval.

- support:

  Optional support constraint for conductance evaluation.

- support_probs:

  Optional length-2 probability vector used with `support = "focal"` to
  compute clamping bounds from focal-cell quantiles.

- clamp_covariates:

  Optional character vector of covariate names to clamp when
  `support = "focal"`.

- ...:

  Additional arguments passed to methods.

## Value

If `x$stack` is non-`NULL` (the default when `saveStack = TRUE` in
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)),
a three-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with layer names:

- `est`:

  Point-estimate conductance at each active cell.

- `lower<Q>`:

  Lower bound of the `quantile*100`% interval (e.g. `lower95` for the
  default `quantile = 0.95`).

- `upper<Q>`:

  Upper bound of the interval.

Otherwise, a numeric matrix with one row per active cell and columns
`est`, `lower<Q>`, and `upper<Q>`.

Conductance is relative: a common multiplier is absorbed by the
measurement model, so compare values among cells of one surface instead
of reading them as absolute rates.

## See also

[`radish_compatibility`](https://wpeterman.github.io/terradish/reference/radish_compatibility.md)
