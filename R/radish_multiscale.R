.validate_multiscale_covariates <- function(covariates, caller = "`terradish_multiscale()`")
{
  covariates <- .as_spatraster(covariates)
  is_factor_layer <- vapply(seq_len(nlyr(covariates)),
                            function(i) is.factor(covariates[[i]]),
                            logical(1))
  if (any(is_factor_layer))
    stop(caller, " currently supports continuous rasters only")
  covariates
}

.aggregate_covariates <- function(covariates, factor, aggregate_fun = mean)
{
  factor <- as.integer(factor)[1]
  if (is.na(factor) || factor < 1L)
    stop("`factor` must be a positive integer")
  if (factor == 1L)
    return(covariates)
  aggregate(covariates, fact = factor, fun = aggregate_fun, na.rm = TRUE)
}

.normalize_coarse_raster_control <- function(control, data = NULL)
{
  control <- if (is.null(control)) list() else as.list(control)
  defaults <- list(
    factor = 2L,
    aggregate_fun = mean,
    directions = if (!is.null(data$directions)) data$directions else 8L,
    exact_refine = TRUE,
    refine_control = NULL
  )
  control <- modifyList(defaults, control)
  control$factor <- sort(unique(as.integer(control$factor)), decreasing = TRUE)
  if (!length(control$factor) || anyNA(control$factor) || any(control$factor < 1L))
    stop("`approximation_control$factor` must contain positive integers")
  if (!is.function(control$aggregate_fun))
    stop("`approximation_control$aggregate_fun` must be a function")
  control$directions <- as.integer(control$directions)[1]
  if (is.na(control$directions) || !control$directions %in% c(4L, 8L))
    stop("`approximation_control$directions` must be 4 or 8")
  if (!is.null(control$refine_control) && !is.list(control$refine_control))
    stop("`approximation_control$refine_control` must be an optimizer control list")
  control
}

.coarse_raster_surface <- function(data, approximation_control = NULL)
{
  control <- .normalize_coarse_raster_control(approximation_control, data = data)
  coarse_factors <- control$factor[control$factor > 1L]
  if (!length(coarse_factors))
  {
    return(list(surface = data,
                stages = list(),
                control = control,
                used = FALSE,
                duplicate_demes = 0L,
                unique_demes = length(unique(data$demes)),
                full_vertices = nrow(data$x),
                coarse_vertices = nrow(data$x)))
  }

  if (is.null(data$stack))
    stop("`approximation = \"coarse_raster\"` requires `data` to retain its raster stack. Recreate it with `conductance_surface(..., saveStack = TRUE)`.",
         call. = FALSE)

  coords <- .deme_coordinates(data)
  if (is.null(coords))
    stop("Could not reconstruct focal coordinates from `data` for coarse-raster screening.",
         call. = FALSE)

  covariates <- .validate_multiscale_covariates(
    data$stack,
    caller = "`approximation = \"coarse_raster\"`"
  )
  stages <- lapply(coarse_factors, function(fact) {
    coarse_covariates <- .aggregate_covariates(covariates,
                                               factor = fact,
                                               aggregate_fun = control$aggregate_fun)
    coarse_surface <- conductance_surface(coarse_covariates,
                                          coords,
                                          directions = control$directions,
                                          saveStack = TRUE)
    list(
      factor = fact,
      surface = coarse_surface,
      duplicate_demes = length(coarse_surface$demes) - length(unique(coarse_surface$demes)),
      unique_demes = length(unique(coarse_surface$demes)),
      coarse_vertices = nrow(coarse_surface$x)
    )
  })
  names(stages) <- paste0("factor_", coarse_factors)

  list(surface = stages[[1L]]$surface,
       stages = stages,
       control = control,
       used = TRUE,
       duplicate_demes = vapply(stages, `[[`, integer(1), "duplicate_demes"),
       unique_demes = vapply(stages, `[[`, integer(1), "unique_demes"),
       full_vertices = nrow(data$x),
       coarse_vertices = vapply(stages, `[[`, integer(1), "coarse_vertices"))
}

#' Coarse-to-fine optimization for large rasters
#'
#' Fits the same resistance-surface model on a sequence of aggregated rasters,
#' using the estimate from each coarser raster as the starting value for the
#' next finer raster.
#'
#' @param formula A formula with a genetic-distance matrix on the left-hand
#'   side and raster covariates on the right-hand side.
#' @param covariates A \code{SpatRaster} of covariates on the finest grid.
#' @param coords Focal coordinates passed to \code{\link{conductance_surface}}.
#' @param factors Integer aggregation factors, ordered from coarse to fine.
#'   \code{1} means the original raster resolution.
#' @param directions Neighborhood definition passed to
#'   \code{\link{conductance_surface}}.
#' @param aggregate_fun Aggregation function for continuous rasters.
#' @param save_surfaces If \code{TRUE}, retain the intermediate
#'   \code{terradish_graph} objects.
#' @param ... Additional arguments passed to \code{\link{terradish}}.
#'
#' @details This helper is aimed at large continuous rasters, where a good
#' coarse-resolution starting point can reduce the amount of full-resolution
#' optimization work. Factor-valued rasters are not currently supported in the
#' aggregation path.
#'
#' For most new analyses, the integrated
#' \code{terradish(..., approximation = "coarse_raster")} interface is usually
#' simpler because it stores approximation metadata directly on the returned
#' model and can optionally refine the final estimate on the full-resolution
#' graph. \code{terradish_multiscale()} remains useful when you want to inspect
#' or retain every intermediate fit explicitly.
#'
#' @return A fitted \code{terradish} object from the finest raster, with an
#' additional \code{$multiscale} component containing the per-level fits.
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#'
#' fit <- terradish_multiscale(
#'   melip.Fst ~ altitude + forestcover,
#'   covariates = covariates,
#'   coords = melip.coords,
#'   factors = c(2, 1),
#'   conductance_model = loglinear_conductance,
#'   measurement_model = generalized_wishart,
#'   nu = 1000,
#'   control = NewtonRaphsonControl(maxit = 3, verbose = FALSE)
#' )
#' fit$multiscale$factors
#'
#' @importFrom terra aggregate
#' @export
terradish_multiscale <- function(formula,
                              covariates,
                              coords,
                              factors = c(4L, 2L, 1L),
                              directions = 8L,
                              aggregate_fun = mean,
                              save_surfaces = FALSE,
                              ...)
{
  stopifnot(inherits(formula, "formula"))
  covariates <- .validate_multiscale_covariates(covariates)
  factors <- sort(unique(as.integer(factors)), decreasing = TRUE)
  if (!length(factors) || anyNA(factors) || any(factors < 1L))
    stop("`factors` must contain positive integers")

  terms_obj <- terms(formula)
  response_idx <- attr(terms_obj, "response")
  if (!response_idx)
    stop("`formula` must have a genetic distance matrix on the left-hand side")
  response_name <- all.vars(formula[[2]])
  if (length(response_name) != 1L)
    stop("The left-hand side of `formula` must be a single matrix object")
  response_value <- get(response_name, parent.frame())
  assign(response_name, response_value, envir = environment())

  theta_start <- NULL
  fits <- vector("list", length(factors))
  surfaces <- if (save_surfaces) vector("list", length(factors)) else NULL

  for (i in seq_along(factors))
  {
    fact <- factors[[i]]
    covariates_i <- .aggregate_covariates(covariates,
                                          factor = fact,
                                          aggregate_fun = aggregate_fun)

    surface_i <- conductance_surface(covariates_i, coords, directions = directions)
    fit_i <- terradish(formula, data = surface_i, theta = theta_start, ...)
    fits[[i]] <- fit_i
    if (save_surfaces)
      surfaces[[i]] <- surface_i
    if (!is.null(fit_i$mle$theta))
      theta_start <- fit_i$mle$theta
  }

  out <- fits[[length(fits)]]
  out$multiscale <- list(
    factors = factors,
    fits = fits,
    surfaces = if (save_surfaces) surfaces else NULL
  )
  out
}

#' Legacy radish multiscale wrapper
#'
#' Deprecated compatibility wrapper retained for older code that still calls
#' \code{radish_multiscale()}.
#'
#' @param ... Additional arguments passed to \code{\link{terradish}}.
#' @name legacy_radish_multiscale_wrapper
#' @keywords internal
NULL

#' @rdname legacy_radish_multiscale_wrapper
#' @export
radish_multiscale <- function(...)
{
  .terradish_deprecate("radish_multiscale", "terradish_multiscale")
  .terradish_forward_call(match.call(), "terradish_multiscale")
}
