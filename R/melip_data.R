#' Melipona Example Data
#'
#' Example data bundled with \pkg{terradish} for demonstrating raster-based
#' conductance modeling. The objects are distributed separately in the package
#' data so they can be loaded individually or together via \code{data(melip)}.
#'
#' @format A collection of objects:
#' \describe{
#'   \item{\code{melip.altitude}}{A wrapped \code{terra::SpatRaster} giving altitude.}
#'   \item{\code{melip.forestcover}}{A wrapped \code{terra::SpatRaster} giving forest cover.}
#'   \item{\code{melip.coords}}{A wrapped \code{terra::SpatVector} of focal-point coordinates.}
#'   \item{\code{melip.Fst}}{A symmetric matrix of pairwise genetic distances among focal points.}
#' }
#'
#' @usage data(melip)
#' @name melip
#' @aliases melip melip.altitude melip.forestcover melip.coords melip.Fst
#' @docType data
#' @keywords datasets
NULL
