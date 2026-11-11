# Crop rasters to a focal-site buffer

Crops a raster stack to the bounding box of a set of focal coordinates
plus a user-supplied map-unit buffer. This is a simple pre-processing
helper for large landscapes where the sampled sites occupy only a small
part of the raster extent.

## Usage

``` r
crop_to_focal_buffer(covariates, coords, buffer)
```

## Arguments

- covariates:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- coords:

  Focal coordinates accepted by
  [`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md):
  a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  matrix, or data frame with coordinate columns.

- buffer:

  Nonnegative buffer distance in the raster's map units. A scalar
  applies the same buffer in both x and y directions; a length-two
  vector is interpreted as `c(x_buffer, y_buffer)`.

## Value

A cropped
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).
Any `terradish_scale` metadata attached by
[`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)
is preserved.

## Details

Cropping reduces the number of graph vertices before any Laplacian
factorization occurs. It is most useful when sampling locations are
spatially clustered inside a much larger raster. It will not help much
when focal sites span the whole raster, because the cropped extent will
be close to the original extent.

The requested buffer is expanded internally by half a cell width in each
direction so the cropped raster safely contains the cells holding the
focal sites even when `buffer = 0`. If the raster is in
longitude/latitude, the buffer is in degrees; use a projected raster for
direct distance interpretation.

[`conductance_surface()`](https://wpeterman.github.io/terradish/reference/conductance_surface.md)
also accepts `crop_buffer = ` and applies this helper before graph
construction. Cropping can increase resistance by removing paths outside
the retained extent. An unpublished development check on one landscape
found increases of 15–23 percent at a two-cell buffer and 1–2 percent at
ten cells. These results are specific to that landscape. Compare fits
over larger buffers as a sensitivity check. Keep the response, focal
sites, covariate scaling, and formula fixed, then refit with
successively larger buffers. Compare coefficients, confidence intervals,
and fitted resistance distances. Report the buffer and whether
increasing it changes the substantive conclusion. A message flags
buffers smaller than half the maximum intersite distance.

## See also

[`conductance_surface`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)

Other graph setup:
[`conductance_surface()`](https://wpeterman.github.io/terradish/reference/conductance_surface.md),
[`scale_covariates()`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)

## Examples

``` r
library(terra)
r <- terra::rast(nrows = 20, ncols = 20, vals = 1:400,
                 xmin = 0, xmax = 20, ymin = 0, ymax = 20)
pts <- matrix(c(4.5, 4.5,
                6.5, 6.5,
                7.5, 5.5),
              ncol = 2, byrow = TRUE)
cropped <- crop_to_focal_buffer(r, pts, buffer = 2)
cropped
#> class       : SpatRaster
#> size        : 7, 8, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 2, 10, 2, 9  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (CRS84) (OGC:CRS84)
#> source(s)   : memory
#> name        : lyr.1
#> min value   :   223
#> max value   :   350

# The same crop can be applied automatically while building the graph.
surface <- conductance_surface(r, pts, directions = 4, crop_buffer = 2)
nrow(surface$x)  # number of retained graph cells
#> [1] 56
```
