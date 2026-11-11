# Assign sampling sites to cross-validation folds

Creates either spatially blocked folds by k-means clustering projected
site coordinates or balanced random folds. Fold labels can also be
supplied directly to
[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)
when study-specific blocks are more defensible than an automated
partition.

## Usage

``` r
terradish_folds(
  pts,
  k = 5L,
  method = c("spatial_kmeans", "random"),
  seed = 1L,
  projected = NULL,
  nstart = 25L
)
```

## Arguments

- pts:

  Focal-site coordinates as a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  matrix, or data frame.

- k:

  Number of folds. Must be at least 2 and smaller than the number of
  sites.

- method:

  Fold construction method. `"spatial_kmeans"` applies k-means to the
  unscaled projected coordinates. `"random"` assigns sites to balanced
  random folds.

- seed:

  Optional random seed.

- projected:

  Logical indicating whether the coordinates are projected. For spatial
  k-means with a `SpatVector`, `NULL` uses its coordinate reference
  system. For spatial k-means with a matrix or data frame, this argument
  must be supplied explicitly. It is ignored for random folds.

- nstart:

  Number of random starts used by k-means.

## Value

An integer vector of fold labels with attributes recording the
construction method, seed, and number of folds.

## Details

Coordinates are deliberately not standardized before spatial clustering.
Consequently, Euclidean separation is measured in the map units of the
projected coordinate reference system. Automated k-means blocks are a
reproducible diagnostic, not a substitute for ecological blocking based
on watersheds, barriers, sampling campaigns, or other design
information.

## See also

[`terradish_cv_folds`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md)

Other inference and model comparison:
[`aic_table()`](https://wpeterman.github.io/terradish/reference/aic_table.md),
[`mlpe_response_change()`](https://wpeterman.github.io/terradish/reference/mlpe_response_change.md),
[`terradish_cv_folds()`](https://wpeterman.github.io/terradish/reference/terradish_cv_folds.md),
[`terradish_ibe_ratio()`](https://wpeterman.github.io/terradish/reference/terradish_ibe_ratio.md),
[`terradish_inference`](https://wpeterman.github.io/terradish/reference/terradish_inference.md),
[`terradish_rescale_nu()`](https://wpeterman.github.io/terradish/reference/terradish_rescale_nu.md)

## Examples

``` r
data(melip)
coords <- terra::unwrap(melip.coords)

# The melip coordinates are longitude and latitude. Project them to a local
# azimuthal equidistant system centered on the sites, so that k-means
# distances are in meters.
lonlat <- terra::crds(coords)
local_crs <- sprintf(
  "+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
  mean(lonlat[, 2]), mean(lonlat[, 1]))
coords_projected <- terra::project(coords, local_crs)

folds <- terradish_folds(coords_projected, k = 3, seed = 1)
table(folds)   # number of sites in each spatial block
#> folds
#>  1  2  3 
#> 18  9 10 
```
