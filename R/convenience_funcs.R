#' Rescale values to the unit interval
#'
#' Rescales numeric values so the minimum non-missing value is 0 and the
#' maximum is 1. \code{terra::SpatRaster} inputs are rescaled in place and
#' returned as rasters.
#'
#' @param x Numeric vector, matrix, dist object, or \code{terra::SpatRaster}.
#'
#' @return An object of the same general type as \code{x}, rescaled to
#'   \code{[0, 1]} where possible.
#'
#' @examples
#' scale_to_0_1(c(2, 4, 6))
#'
#' library(terra)
#' r <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 3, 5, NA))
#' scale_to_0_1(r)
#'
#' @export
scale_to_0_1 <- function(x)
{
  if (inherits(x, "PackedSpatRaster"))
    x <- unwrap(x)

  if (inherits(x, "SpatRaster"))
  {
    vals <- values(x, dataframe = FALSE)
    rng <- range(vals, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0)
    {
      vals[!is.na(vals)] <- 0
    }
    else
    {
      vals <- (vals - rng[1]) / diff(rng)
    }
    values(x) <- vals
    return(x)
  }

  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0)
  {
    out <- x
    out[!is.na(out)] <- 0
    return(out)
  }

  (x - rng[1]) / diff(rng)
}

#' Extract the lower triangle of a matrix
#'
#' Returns the lower-triangular entries of a square matrix as a vector.
#'
#' @param x A square matrix.
#'
#' @return A vector containing the strict lower triangle of \code{x}.
#'
#' @examples
#' lower(matrix(1:9, nrow = 3))
#'
#' @export
lower <- function(x)
{
  stopifnot(is.matrix(x), nrow(x) == ncol(x))
  x[lower.tri(x)]
}

#' Pairwise genetic distances from PCA scores
#'
#' Computes Euclidean distances among individuals using retained principal
#' component axes from an \code{adegenet::genind} object.
#'
#' @param gi A \code{genind} object from \pkg{adegenet}.
#' @param n_axes Number of principal component axes to retain.
#' @param scale Should the resulting distance matrix be rescaled to \code{[0,1]}?
#'
#' @return A symmetric distance matrix.
#'
#' @examples
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   data(nancycats, package = "adegenet")
#'   pca_dist(nancycats, n_axes = 4)
#' }
#'
#' @export
pca_dist <- function(gi, n_axes = 64, scale = TRUE)
{
  a_tab <- .adegenet_tab(gi, NA.method = "mean")
  pc <- prcomp(a_tab)
  n_axes <- min(as.integer(n_axes), ncol(pc$x))
  if (n_axes < 1)
    stop("`n_axes` must retain at least one principal component axis")

  pc_dist <- as.matrix(dist(pc$x[, seq_len(n_axes), drop = FALSE]))
  if (isTRUE(scale))
    pc_dist <- as.matrix(scale_to_0_1(pc_dist))

  pc_dist
}
