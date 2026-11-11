# Build the raster graph from covariates and sampling sites

Reads spatial covariates and focal-point coordinates and builds a
weighted graph suitable for computing conductance-weighted resistance
distances. This is the first step of every `terradish` workflow.

## Usage

``` r
conductance_surface(
  covariates,
  coords,
  directions = 8,
  saveStack = TRUE,
  crop_buffer = NULL
)
```

## Arguments

- covariates:

  A `SpatRaster` (or `PackedSpatRaster`) whose layers contain the
  spatial covariates to be used in the conductance model. Layer names
  must match the right-hand side of the formula passed to
  [`terradish`](https://wpeterman.github.io/terradish/reference/terradish.md).
  Name the layers first, then scale continuous layers with
  [`scale_covariates`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)
  before calling this function.

- coords:

  Locations of the focal sampling sites, in the same coordinate
  reference system as `covariates`. Supported inputs:
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  of points, a two-column numeric matrix (x in column 1, y in column 2),
  or a data frame with columns named `x` and `y` (or at least two
  columns, using the first two).

- directions:

  Adjacency rule for the graph. `4` connects each cell to its four
  horizontal/vertical neighbors (rook adjacency); `8` additionally
  includes the four diagonal neighbors (queen adjacency). **`8` is
  recommended** for landscape genetics: it allows diagonal movement,
  produces smoother resistance surfaces, and is less sensitive to grid
  orientation artifacts.

- saveStack:

  Logical. If `TRUE` (default), the masked `SpatRaster` is stored inside
  the returned `terradish_graph` object. Required for raster-based
  conductance visualization (`plot(fit, type = "surface")`) and
  marginal-association plots. Set to `FALSE` to reduce memory use when
  the raster is large and visualization is not needed.

- crop_buffer:

  Optional nonnegative map-unit buffer to crop `covariates` around the
  focal coordinates before graph construction. When supplied,
  [`crop_to_focal_buffer`](https://wpeterman.github.io/terradish/reference/crop_to_focal_buffer.md)
  is called internally. A scalar applies the same buffer in both x and y
  directions; a length-two vector gives `c(x_buffer, y_buffer)`. Units
  are the raster's map units (degrees for lon/lat, meters for projected
  rasters). See Details.

## Value

An object of class `"terradish_graph"` (also inheriting
`"radish_graph"`) containing:

- `demes`:

  Integer indices of cells corresponding to focal sites.

- `x`:

  Data frame of covariate values at non-missing cells.

- `adj`:

  Edge list (upper-triangular, 0-based).

- `edge_pairs`:

  Integer matrix of adjacent cell pairs.

- `vertex_coordinates`:

  xy-coordinate matrix of active cells.

- `directions`:

  The adjacency rule (4 or 8).

- `covariates`:

  Names of the covariate columns.

- `laplacian`, `rhs`, `reduced_index`, `choleski_templates`, `choleski`:

  Internal solver state: the reduced graph Laplacian template, the
  right-hand sides for the focal sites, the retained vertex index, and
  symbolic sparse Cholesky factorizations. These are reused by the
  fitting functions and are not meant to be modified.

- `stack`:

  The masked `SpatRaster` (or `NULL` if `saveStack = FALSE`).

## Details

**Missing values:** `NA` is propagated uniformly across all layers so
every cell is either complete or entirely missing. A warning is issued
if some cells are `NA` in only a subset of layers, because this
indicates that different covariates cover different areas; model
comparison across raster stacks with different `NA` patterns is
problematic since the graph vertices differ. Always compare models using
the same `terradish_graph` object.

**Disconnected components:** after masking, the graph may contain cells
that are isolated from the main component. The function retains only the
largest connected component and warns about the number of pruned cells.
The function aborts if any focal site belongs to a disconnected
component.

**Cropping:** when `crop_buffer` is supplied, cropping changes the graph
domain. Always use the same `terradish_graph` for all models you intend
to compare. Cropping removes possible movement paths and can inflate
resistance. In an unpublished development check on one landscape,
resistance increased by 15–23 percent with a two-cell buffer and by 1–2
percent with ten cells. These are case-specific results; refit with
several buffer values before interpreting a cropped landscape.

**Categorical covariates:** raster layers representing categorical
variables must be stored as factor-valued `SpatRaster` layers (see
[`as.factor`](https://rspatial.github.io/terra/reference/is.bool.html)
and the example below). Level names are taken from the `VALUE` column of
the levels table when present; otherwise the first non-`ID` column is
used. Factor levels are dummy-coded by
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) during model
fitting, with one level as the reference.

**Directions:** eight-neighbor adjacency is the default and is preferred
unless there is a specific reason to restrict movement to the cardinal
directions. `directions = 4` removes the diagonal edges, so it can
overstate resistance (understate connectivity) along diagonal corridors.
Each edge has conductance \\c_i + c_j\\; diagonal edges receive no
additional distance weighting.

**Sites sharing a cell:** two or more focal sites that fall in the same
raster cell trigger a warning, because their resistance distance is
zero. Use a finer raster, or keep one site per cell as in the example.

## See also

Other graph setup:
[`crop_to_focal_buffer()`](https://wpeterman.github.io/terradish/reference/crop_to_focal_buffer.md),
[`scale_covariates()`](https://wpeterman.github.io/terradish/reference/scale_covariates.md)

## Examples

``` r
library(terra)

data(melip)
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
coords <- terra::unwrap(melip.coords)

# Name the layers first, then scale them.
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

# A few melip sites fall in the same raster cell. Keep one site per cell
# (and the matching rows and columns of the genetic matrix).
keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
coords <- coords[keep]
Fst <- melip.Fst[keep, keep]

surface <- conductance_surface(covariates, coords)  # eight neighbors
surface$directions
#> [1] 8
nrow(surface$x)       # number of graph vertices (non-missing cells)
#> [1] 22443
length(surface$demes) # number of focal sites
#> [1] 34

# For large rasters where sites occupy only part of the map, crop first.
surface_cropped <- conductance_surface(covariates, coords,
                                       crop_buffer = 0.05)
#> The crop buffer is less than half the maximum distance between focal sites. Cropping can inflate resistance; compare fits with a larger buffer.

# Categorical covariates must be factor-valued layers (see Details).
forestcover <- terra::unwrap(melip.forestcover)
forestcover_class <- cut(terra::values(forestcover)[, 1],
                         breaks = c(0, 1/6, 1/3, 1))
forestcover_cat <- terra::setValues(forestcover,
                                    as.numeric(forestcover_class))
forestcover_cat <- terra::as.factor(forestcover_cat)
RAT <- levels(forestcover_cat)[[1]]
RAT$VALUE <- levels(forestcover_class) # explicitly define level names
levels(forestcover_cat) <- RAT
names(forestcover_cat) <- "forestcover"

surface_cat <- conductance_surface(forestcover_cat, coords)
levels(surface_cat$x$forestcover)
#> NULL
```
