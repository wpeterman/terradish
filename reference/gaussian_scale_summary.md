# Summarize fitted Gaussian smoothing scales

Converts fitted Gaussian `sigma` parameters into map-unit, cell-based,
and Gaussian-kernel interpretation summaries for a fitted model produced
with
[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md).

## Usage

``` r
gaussian_scale_summary(
  object,
  probabilities = c(0.5, 0.95),
  distance_per_map_unit = NULL,
  distance_unit = NULL
)
```

## Arguments

- object:

  A fitted
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md)
  object using the Gaussian scale-aware conductance model.

- probabilities:

  Probabilities used to summarize the one-dimensional half-width and
  two-dimensional radial extent of the Gaussian kernel.

- distance_per_map_unit:

  Optional scalar conversion factor used to express `sigma` and the
  derived radii in user-supplied distance units. For example, if the
  raster CRS is in meters, use `distance_per_map_unit = 0.001` and
  `distance_unit = "km"` to report kilometers.

- distance_unit:

  Optional label used when `distance_per_map_unit` is supplied.

## Value

A data frame with one row per fitted `sigma` parameter. `covariate`
identifies the layer. `sigma` is its fitted map-unit scale;
`sigma_lower` and `sigma_upper` are fitting bounds. `sigma_internal`,
`sigma_conversion`, and `sigma_conversion_mode` describe the
optimization units. `native_unit` labels map units or degrees, and
`sigma_cells_x`, `sigma_cells_y` express the scale in cell widths.
Probability-specific `axis_*` and `radial_*` columns give the
theoretical Gaussian distances described above, in native units and,
when requested, converted units. They describe the untruncated Gaussian;
retained raster windows can omit some of that mass.

## Details

The returned table always reports `sigma` on the raster's native map
scale and in raster cell widths. When sigma conversion is enabled during
fitting, the table also reports the internal optimization scale and its
conversion factor back to map units. For each probability `p`, `axis_*`
gives the one-dimensional half-width \\\sigma \Phi^{-1}((1+p)/2)\\,
while `radial_*` gives the isotropic two-dimensional radius
\\\sigma\sqrt{-2\log(1-p)}\\ containing proportion `p` of the Gaussian
kernel mass. Here \\\Phi^{-1}\\ is the standard normal quantile function
and \\\sigma\\ is the smoothing scale. These quantities summarize the
fitted raster-smoothing kernel. They are not estimates of dispersal
distance, movement distance, or home-range size. Use
`summary(object)$sigma_table` to inspect standard errors and bound
flags, and `confint(object)` for bound-truncated intervals; smoothing
scales have no reported z or p values. Use
[`gaussian_scale_profile`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md)
when a profile interval is needed. A positive scale can reflect
dispersal averaging, and even the half-cell lower bound applies some
smoothing.

If the retained raster is in longitude/latitude, the native-unit results
are in degrees. In that case, use a projected raster for direct distance
interpretation, or provide an approximate `distance_per_map_unit`
conversion if you need a quick descriptive summary. Without a conversion
the function warns once per session, and the `native_unit` column reads
`"degrees"` on every call.

## See also

[`gaussian_smoothed_loglinear_conductance`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`gaussian_scale_profile`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md),
[`plot.terradish`](https://wpeterman.github.io/terradish/reference/plot.terradish.md)

Other spatial scale:
[`gaussian_scale_profile()`](https://wpeterman.github.io/terradish/reference/gaussian_scale_profile.md),
[`gaussian_smoothed_loglinear_conductance()`](https://wpeterman.github.io/terradish/reference/gaussian_smoothed_loglinear_conductance.md),
[`terradish_scale_optim()`](https://wpeterman.github.io/terradish/reference/terradish_scale_optim.md)

## Examples

``` r
# \donttest{
data(melip)
forestcover <- terra::unwrap(melip.forestcover)
coords <- terra::unwrap(melip.coords)
names(forestcover) <- "forestcover"

# Coarsen the raster so the example runs quickly; the Gaussian model
# re-smooths the layer at every optimizer step, so cost scales with cell
# count. Use the full resolution in a real analysis.
forestcover <- terra::aggregate(forestcover, fact = 3, na.rm = TRUE)

# Keep one site per raster cell (some melip sites share a cell).
keep <- !duplicated(terra::cellFromXY(forestcover, terra::crds(coords)))
Fst <- melip.Fst[keep, keep]

# Feed the model the raw (unsmoothed, unstandardized) raster: it smooths and
# standardizes internally at each candidate sigma.
surface <- conductance_surface(forestcover, coords[keep], saveStack = TRUE)
fit <- terradish(Fst ~ forestcover,
                 data = surface,
                 conductance_model =
                   gaussian_smoothed_loglinear_conductance(surface),
                 measurement_model = leastsquares,
                 optimizer = "bfgs")

# sigma is in map units. The melip rasters are unprojected, so those units
# are degrees; distance_per_map_unit rescales the reported distances
# (roughly 111 km per degree near the equator). Without a conversion the
# function warns that the scale is in degrees.
scales <- gaussian_scale_summary(fit, distance_per_map_unit = 111,
                                 distance_unit = "km")
scales
#>     covariate     sigma sigma_lower sigma_upper sigma_internal sigma_conversion
#> 1 forestcover 0.2621598      0.0675        1.17       1.941924            0.135
#>   sigma_conversion_mode native_unit sigma_cells_x sigma_cells_y   axis_50
#> 1                  cell     degrees      1.941924      1.941924 0.1768241
#>   radial_50   axis_95 radial_95 sigma_km sigma_lower_km sigma_upper_km
#> 1 0.3086696 0.5138238 0.6417008 29.09974         7.4925         129.87
#>   axis_50_km radial_50_km axis_95_km radial_95_km
#> 1   19.62748     34.26232   57.03444     71.22879
# One row per smoothed layer: sigma, its bounds, and the kernel radii.
scales[, c("covariate", "sigma", "sigma_lower", "sigma_upper")]
#>     covariate     sigma sigma_lower sigma_upper
#> 1 forestcover 0.2621598      0.0675        1.17
# }
```
