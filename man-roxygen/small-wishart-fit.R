#' @examples
#' # A small model-based example, with fixed information supplied explicitly.
#' r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
#' terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
#' names(r) <- "x"
#' sites <- terra::xyFromCell(r, c(1, 8, 20, 37, 57, 64))
#' graph <- conductance_surface(r, sites)
#' simulated <- simulate_covariance_response(c(x = 0.6), ~x, graph,
#'                                           nu = 1000, seed = 93)
#' response <- simulated$covariance
#' fit <- terradish(response ~ x, graph, measurement_model = wishart_covariance,
#'                  nu = 1000)
