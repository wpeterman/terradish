# Rescale values to the unit interval

Rescales numeric values so the minimum non-missing value maps to 0 and
the maximum maps to 1.
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
inputs are rescaled layer-by-layer using only the finite (non-`NA`)
values.

## Usage

``` r
scale_to_0_1(x)
```

## Arguments

- x:

  Numeric vector, matrix, `dist` object, or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).
  `PackedSpatRaster` objects (from
  [`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html))
  are unwrapped automatically.

## Value

An object of the same class as `x` with values in `[0, 1]` (non-`NA`
entries only).

## Details

The transformation applied to each element is \\(x - \min) / (\max -
\min)\\, computed over non-`NA` values. `NA` values are propagated
unchanged.

If all non-missing values are identical (constant layer or vector), the
output is set to 0 for all non-missing entries rather than producing
`NaN` from a zero-range denominator.

For rasters, each layer is rescaled independently using its own minimum
and maximum.

## See also

[`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)
for standardized (z-score or min-max) scaling with metadata retained for
back-transformation.

## Examples

``` r
scale_to_0_1(c(2, 4, 6))
#> [1] 0.0 0.5 1.0
scale_to_0_1(c(5, 5, NA))  # constant non-NA values become 0
#> [1]  0  0 NA

library(terra)
r <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 3, 5, NA))
scale_to_0_1(r)
#> class       : SpatRaster
#> size        : 2, 2, 1  (nrow, ncol, nlyr)
#> resolution  : 180, 90  (x, y)
#> extent      : -180, 180, -90, 90  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
#> source(s)   : memory
#> name        : lyr.1
#> min value   :     0
#> max value   :     1
```
