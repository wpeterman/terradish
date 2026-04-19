#' Crop rasters to a focal-site buffer
#'
#' Crops a raster stack to the bounding box of a set of focal coordinates plus
#' a user-supplied map-unit buffer. This is a simple pre-processing helper for
#' large landscapes where the sampled sites occupy only a small part of the
#' raster extent.
#'
#' @param covariates A \code{terra::SpatRaster}.
#' @param coords Focal coordinates accepted by \code{\link{conductance_surface}}:
#'   a \code{terra::SpatVector}, matrix, or data frame with coordinate columns.
#' @param buffer Nonnegative buffer distance in the raster's map units. A scalar
#'   applies the same buffer in both x and y directions; a length-two vector is
#'   interpreted as \code{c(x_buffer, y_buffer)}.
#'
#' @details
#' Cropping reduces the number of graph vertices before any Laplacian
#' factorization occurs. It is most useful when sampling locations are spatially
#' clustered inside a much larger raster. It will not help much when focal
#' sites span the whole raster, because the cropped extent will be close to the
#' original extent.
#'
#' The requested buffer is expanded internally by half a cell width in each
#' direction so the cropped raster safely contains the cells holding the focal
#' sites even when \code{buffer = 0}. If the raster is in longitude/latitude,
#' the buffer is in degrees; use a projected raster for direct distance
#' interpretation.
#'
#' \code{conductance_surface()} also accepts \code{crop_buffer = } and applies
#' this helper before graph construction.
#'
#' @return A cropped \code{terra::SpatRaster}. Any \code{terradish_scale}
#'   metadata attached by \code{\link{scale_covariates}} is preserved.
#'
#' @examples
#' library(terra)
#' r <- terra::rast(nrows = 20, ncols = 20, vals = 1:400)
#' pts <- matrix(c(4.5, 4.5, 6.5, 6.5), ncol = 2, byrow = TRUE)
#' cropped <- crop_to_focal_buffer(r, pts, buffer = 2)
#' cropped
#'
#' # The same crop can be applied automatically while building the graph.
#' surface <- conductance_surface(r, pts, crop_buffer = 2)
#' surface$dim
#'
#' @export
crop_to_focal_buffer <- function(covariates, coords, buffer)
{
  covariates <- .as_spatraster(covariates)
  coords <- .coords_matrix(coords, covariates)

  buffer <- as.numeric(buffer)
  if (!length(buffer) || length(buffer) > 2L ||
      any(!is.finite(buffer)) || any(buffer < 0))
    stop("`buffer` must be one or two finite nonnegative numbers",
         call. = FALSE)
  if (length(buffer) == 1L)
    buffer <- rep(buffer, 2L)

  resolution <- abs(res(covariates[[1]]))
  pad_x <- buffer[[1]] + resolution[[1]] / 2
  pad_y <- buffer[[2]] + resolution[[2]] / 2

  crop_extent <- ext(
    min(coords[, 1]) - pad_x,
    max(coords[, 1]) + pad_x,
    min(coords[, 2]) - pad_y,
    max(coords[, 2]) + pad_y
  )

  scale_meta <- attr(covariates, "terradish_scale", exact = TRUE)
  out <- crop(covariates, crop_extent)
  if (!is.null(scale_meta))
    attr(out, "terradish_scale") <- scale_meta
  out
}
