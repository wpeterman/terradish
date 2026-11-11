# Scale raster covariates and retain original-scale metadata

Transforms each layer of a `SpatRaster` and stores the original-scale
offset and divisor as a `"terradish_scale"` attribute. When the scaled
raster is used to build a
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
marginal-association plots
(`plot(fit, type = "marginal", data = surface)`) can automatically
back-transform x-axes to the original covariate units.

## Usage

``` r
scale_covariates(
  covariates,
  method = c("zscore", "minmax"),
  center = TRUE,
  scale = TRUE,
  reference = NULL
)
```

## Arguments

- covariates:

  A `SpatRaster` or `PackedSpatRaster` whose layers should be scaled.

- method:

  Scaling method applied to each layer independently. `"zscore"`
  (default) subtracts the layer mean and divides by the standard
  deviation, giving zero mean and unit variance. `"minmax"` maps the
  minimum to 0 and the maximum to 1.

- center:

  Logical. Only used when `method = "zscore"`. If `TRUE` (default),
  subtract the layer mean before dividing by the standard deviation.

- scale:

  Logical. Only used when `method = "zscore"`. If `TRUE` (default),
  divide by the layer standard deviation. If `FALSE` with
  `center = TRUE`, layers are centered but not divided.

- reference:

  A previous scaled raster or its `terradish_scale` attribute. Applies
  the stored centers and scales by layer name instead of estimating new
  values. When supplied, `method`, `center`, and `scale` do not change
  the stored transformation.

## Value

A `SpatRaster` with transformed values and a `"terradish_scale"`
attribute.

## Details

Scaling is strongly recommended before fitting conductance models. In
the loglinear model, conductance is \\\exp(\theta \cdot x)\\. When `x`
has a large range (e.g. altitude in meters), a small change in
\\\theta\\ produces enormous conductance differences, making the
objective surface poorly conditioned and the optimizer slow to converge.

Missing values (`NA`) are excluded from the computation of the mean,
standard deviation, minimum, and range, and are left as `NA` in the
output.

If a layer is constant (all values identical), the scale is set to 1 to
avoid division by zero; the layer is then returned unchanged (or
centered to 0 if `center = TRUE`).

The `"terradish_scale"` attribute is a named list with one element per
layer; each element is a named numeric vector with entries `center` and
`scale` (the subtracted offset and the divisor, respectively).
`crop_to_focal_buffer` and `conductance_surface` preserve this
attribute.

Name the layers before scaling. The stored centers and scales are keyed
by layer name, so renaming a raster after scaling leaves stale names in
the `"terradish_scale"` attribute and breaks back-transformation.

## See also

[`scale_to_0_1`](https://wpeterman.github.io/terradish/reference/scale_to_0_1.md),
[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)

Other graph setup:
[`conductance_surface()`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`crop_to_focal_buffer()`](https://wpeterman.github.io/terradish/reference/crop_to_focal_buffer.md)

## Examples

``` r
library(terra)
r <- c(terra::rast(nrows = 2, ncols = 2, vals = c(100, 200, 300, NA)),
       terra::rast(nrows = 2, ncols = 2, vals = c(0.2, 0.5, 0.4, NA)))
names(r) <- c("altitude", "forestcover")

scaled_r <- scale_covariates(r)
attr(scaled_r, "terradish_scale")  # per-layer center and scale
#> $altitude
#> center  scale 
#>    200    100 
#> 
#> $forestcover
#>    center     scale 
#> 0.3666667 0.1527525 
#> 

# min-max scaling to [0, 1]
scaled_mm <- scale_covariates(r, method = "minmax")

# Preserve the training transformation on new values, matched by layer name.
new_scaled <- scale_covariates(r + 10, reference = scaled_r)
terra::values(new_scaled)
#>      altitude forestcover
#> [1,]     -0.9    64.37428
#> [2,]      0.1    66.33824
#> [3,]      1.1    65.68358
#> [4,]       NA          NA
# Values can fall outside the training range; they are not re-centered.
```
