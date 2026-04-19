.as_spatraster <- function(covariates)
{
  if (inherits(covariates, "SpatRaster"))
    return(covariates)
  if (inherits(covariates, "PackedSpatRaster"))
    return(unwrap(covariates))
  if (inherits(covariates, c("Raster", "RasterStack", "RasterBrick", "RasterLayer")))
    return(rast(covariates))
  stop("'covariates' must be a terra SpatRaster")
}

.coords_matrix <- function(coords, template)
{
  if (inherits(coords, "PackedSpatVector"))
    coords <- unwrap(coords)
  if (inherits(coords, "SpatVector"))
    return(geom(coords)[, c("x", "y"), drop = FALSE])
  if (is.matrix(coords))
  {
    if (ncol(coords) < 2)
      stop("'coords' matrix must have at least two columns")
    return(coords[, 1:2, drop = FALSE])
  }
  if (is.data.frame(coords))
  {
    if (all(c("x", "y") %in% names(coords)))
      return(as.matrix(coords[, c("x", "y"), drop = FALSE]))
    if (ncol(coords) < 2)
      stop("'coords' data.frame must have x/y columns or at least two columns")
    return(as.matrix(coords[, 1:2, drop = FALSE]))
  }
  stop("'coords' must be a SpatVector, matrix, or data.frame")
}

.factor_levels <- function(x)
{
  lev <- levels(x)[[1]]
  ids <- lev$ID
  label_col <- if ("VALUE" %in% names(lev)) "VALUE" else setdiff(names(lev), "ID")[1]
  list(ids = ids, labels = lev[[label_col]])
}

.reduced_rhs <- function(n_vertices, demes)
{
  stopifnot(n_vertices >= 2L)
  rhs <- matrix(-1 / n_vertices, nrow = n_vertices - 1L, ncol = length(demes))
  keep <- demes < n_vertices
  if (any(keep))
  {
    idx <- cbind(demes[keep], which(keep))
    rhs[idx] <- rhs[idx] + 1
  }
  rhs
}

#' Create a parameterized conductance surface
#'
#' Given a set of spatial covariates and a set of spatial coordinates, create a
#' graph representing a parameterized conductance surface.
#'
#' @param covariates A \code{SpatRaster} containing spatial covariates
#' @param coords A point object containing coordinates for a set of focal cells,
#'   with the same projection as \code{covariates}. Supported inputs are
#'   \code{SpatVector}, matrices, and data frames with \code{x}/\code{y}
#'   columns.
#' @param directions If \code{4}, consider only horizontal/vertical neighbours as adjacent; if \code{8}, also consider diagonal neighbours as adjacent
#' @param saveStack If \code{TRUE}, the \code{SpatRaster} is returned with missing data masked uniformly across layers
#' @param crop_buffer Optional nonnegative map-unit buffer around the focal
#'   coordinates. When supplied, \code{covariates} are cropped with
#'   \code{\link{crop_to_focal_buffer}} before the graph is constructed. This
#'   can substantially reduce graph size for large rasters where focal sites
#'   occupy only part of the landscape. A scalar uses the same buffer in the
#'   x and y directions; a length-two vector is interpreted as
#'   \code{c(x_buffer, y_buffer)}.
#'
#' @details NAs are shared across raster layers in \code{covariates}, and a
#' warning is thrown if a given cell has mixed NA and non-NA values across the
#' stack. Comparing models with different patterns of missing spatial data
#' (e.g. fit to different stacks of rasters) can give superficially
#' inconsistant results, as these essentially involve different sets of
#' vertices. Thus model comparison should use models fitted to the same
#' \code{terradish_graph} object.
#'
#' Disconnected components are identified and removed, so that only the largest connected component in the graph is retained. The function aborts if there are focal cells that belong to a disconnected component.
#'
#' If \code{crop_buffer} is supplied, cropping happens before missing-cell
#' masking and graph construction. This changes the graph domain, so use the
#' same cropped \code{terradish_graph} object for all models you intend to
#' compare. A buffer that is too small can omit landscape context that may
#' affect resistance distances; for sensitivity analyses, refit with several
#' buffer widths and confirm that coefficients and likelihoods are stable.
#'
#' Categorical raster layers should be stored as factor-valued
#' \code{SpatRaster} layers, see \code{\link[terra]{as.factor}} and the
#' examples below. The names of levels are taken from the \code{VALUE} column
#' of the associated levels table when present; otherwise the first non-\code{ID}
#' column is used.
#'
#' @seealso \code{\link{terradish}}, \code{\link[terra]{rast}}
#'
#' @return An object of class \code{terradish_graph}
#'
#' @references
#'
#' Pope NS. In prep. Fast gradient-based optimization of resistance surfaces.
#'
#' @examples
#'
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
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' # For large rasters where sites occupy only part of the map, crop first.
#' surface_cropped <- conductance_surface(
#'   covariates,
#'   melip.coords,
#'   directions = 8,
#'   crop_buffer = 0.05
#' )
#'
#' # categorical covariates:
#' # raster layers of categorical covariates should be factor-valued, see 'details'
#' forestcover_class <- cut(terra::values(melip.forestcover)[,1], breaks = c(0, 1/6, 1/3, 1))
#' melip.forestcover_cat <- terra::setValues(melip.forestcover, as.numeric(forestcover_class))
#' melip.forestcover_cat <- terra::as.factor(melip.forestcover_cat)
#'
#' RAT <- levels(melip.forestcover_cat)[[1]]
#' RAT$VALUE <- levels(forestcover_class) # explicitly define level names
#' levels(melip.forestcover_cat) <- RAT
#'
#' covariates_cat <- c(melip.forestcover_cat, melip.altitude)
#' names(covariates_cat) <- c("forestcover", "altitude")
#'
#' @export

conductance_surface <- function(covariates, coords, directions=4, saveStack=TRUE,
                                crop_buffer = NULL)
{
  covariates <- .as_spatraster(covariates)
  stopifnot(directions %in% c(4, 8))

  coords <- .coords_matrix(coords, covariates)
  if (!is.null(crop_buffer))
    covariates <- crop_to_focal_buffer(covariates, coords, buffer = crop_buffer)
  if (nlyr(covariates) < 1)
    stop("'covariates' must contain at least one layer")

  # throw a warning if missing cells are not identical across layers
  spdat <- values(covariates, dataframe = FALSE)
  missing_count <- rowSums(is.na(spdat))
  missing <- missing_count > 0
  if (any(missing_count > 0 & missing_count < ncol(spdat)))
    warning("Missing cells are not identical across rasters; be careful regarding model selection (see ?conductance_surface)")

  # share missing cells across layers and ensure graph is fully connected
  if (any(missing))
  {
    spdat[missing, ] <- NA_real_
    cr <- covariates[[1]]
    values(cr) <- ifelse(missing, NA_real_, 1)
    cr <- patches(cr, directions = directions)
    components <- values(cr, dataframe = FALSE)[,1]
    connected_component <- as.integer(names(which.max(table(components))))
    disconnected <- components != connected_component
    disconnected[is.na(disconnected)] <- FALSE
    if (any(disconnected))
      warning(paste("Pruned", sum(disconnected), "disconnected cells across rasters"))
    missing <- disconnected | missing
    spdat[missing, ] <- NA_real_
    values(covariates) <- spdat
  }

  # get adjacency list
  active_cells <- which(!missing)
  adj <- adjacent(covariates[[1]], cells = active_cells, directions = directions,
                  pairs = TRUE, include = FALSE)
  adj <- adj[adj[,1] < adj[,2], , drop = FALSE]

  # find cells of demes
  cells <- unmapped_cells <- cellFromXY(covariates[[1]], coords)
  if (any(is.na(unmapped_cells)))
    stop("At least one deme is located outside the raster extent")

  # remove NAs and remap indices of adjacency list to be contiguous
  cell_map <- integer(ncell(covariates[[1]]))
  cell_map[active_cells] <- seq_along(active_cells)
  vertex_coordinates <- xyFromCell(covariates[[1]], active_cells)
  spdat <- spdat[active_cells, , drop = FALSE]
  adj <- cbind(cell_map[adj[,1]], cell_map[adj[,2]])
  adj <- adj[adj[,1] > 0 & adj[,2] > 0, , drop = FALSE]
  cells <- cell_map[cells]

  # check that cells lie on connected portion of raster
  if (any(missing[unmapped_cells]))
    stop("At least one deme is located on a missing cell")

  # figure out which raster layers are factors
  is_factor <- vapply(seq_len(nlyr(covariates)),
    function(i) isTRUE(is.factor(covariates[[i]])),
                      logical(1))
  spdat <- as.data.frame(spdat)
  if (any(is_factor))
  {
    factors <- names(covariates)[is_factor]
    warning("Treating covariates \"", paste(factors, collapse = "\" \""), "\" as factors")
    for (i in factors)
    {
      lev <- .factor_levels(covariates[[i]])
      spdat[, i] <- factor(lev$labels[match(spdat[, i], lev$ids)])
    }
  }

  # form and factorize Laplacian
  N    <- nrow(spdat)
  edge_pairs <- adj
  storage.mode(edge_pairs) <- "integer"
  Q    <- sparseMatrix(i = adj[,1], j = adj[,2], dims = c(N, N),
                       x = -rep(1, nrow(adj)),
                       use.last.ij = TRUE)
  Q    <- forceSymmetric(Q)
  Qd   <- Diagonal(N, x = -rowSums(Q))
  In   <- Diagonal(N)[-N,]
  Qn   <- forceSymmetric(In %*% (Q + Qd) %*% t(In))
  adj  <- rbind(Q@i, rep(1:Q@Dim[2] - 1, diff(Q@p))) #upper-triangular, 0-based
  LQn  <- Cholesky(Qn, LDL = TRUE) 
  rhs  <- .reduced_rhs(N, cells)

  out <- list("demes"        = cells,
              "x"            = spdat,
              "adj"          = adj,
              "edge_pairs"   = edge_pairs,
              "vertex_coordinates" = vertex_coordinates,
              "directions"   = directions,
              "covariates"   = colnames(spdat), 
              "laplacian"    = Q,
              "rhs"          = rhs,
              "reduced_index" = seq_len(N - 1L),
              "choleski_templates" = list(simplicial_ldl = LQn),
              "choleski"     = LQn,
              "stack"        = if(saveStack) covariates else NULL)
  class(out) <- c("terradish_graph", "radish_graph")
  out
}

#' Estimate conductance from a fitted model
#'
#' Returns fitted conductance values, and confidence intervals when available,
#' for a \code{terradish_graph} and fitted \code{terradish} object.
#'
#' @param x A \code{terradish_graph} object.
#' @param fit A fitted object returned by \code{\link{terradish}}.
#' @param quantile Confidence level used to compute conductance intervals.
#' @param ... Additional arguments passed to methods.
#'
#' @return If the graph retains its raster stack, a three-layer
#'   \code{terra::SpatRaster} with estimate and interval bounds. Otherwise a
#'   matrix with one row per active cell.
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
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' fit <- terradish(melip.Fst ~ forestcover + altitude, surface,
#'               loglinear_conductance, leastsquares,
#'               control = NewtonRaphsonControl(maxit = 5))
#'
#' cond <- conductance(surface, fit)
#' cond
#'
#' @export
conductance <- function(x, ...)
{
  UseMethod("conductance")
}

.conductance_graph_impl <- function(x, fit, quantile = 0.95, ...)
{
  stopifnot(inherits(fit, c("terradish", "radish")))
  stopifnot(!fit$fit$boundary && !is.null(fit$mle$theta))

  conductance <- fit$submodels$f(fit$mle$theta)
  ci <- conductance$confint(theta = fit$mle$theta, vcov = -solve(fit$mle$hessian), 
                            quantile = quantile, scale = "conductance")
  colnames(ci) <- paste0(c("lower", "upper"), round(100*quantile, 1))

  if (!is.null(x$stack))
  {
    template <- x$stack[[1]]
    template_values <- values(template, dataframe = FALSE)[,1]
    missing <- is.na(template_values)
    template_values[!missing] <- 1
    values(template) <- template_values
    lower <- upper <- est <- template
    est_values <- lower_values <- upper_values <- template_values
    est_values[!missing] <- conductance$conductance
    lower_values[!missing] <- ci[,1]
    upper_values[!missing] <- ci[,2]
    values(est) <- est_values
    values(lower) <- lower_values
    values(upper) <- upper_values
    out <- c(est, lower, upper)
    names(out) <- c("est", colnames(ci))
    out
  }
  else
  {
    warning("No rasters associated with graph, returning conductance as vector")
    out <- cbind("est" = conductance$conductance, ci)
    out
  }
}

#' @rdname conductance
#' @method conductance terradish_graph
#' @export
conductance.terradish_graph <- function(x, fit, quantile = 0.95, ...)
{
  .conductance_graph_impl(x = x, fit = fit, quantile = quantile, ...)
}

#' @rdname conductance
#' @method conductance radish_graph
#' @export
conductance.radish_graph <- function(x, fit, quantile = 0.95, ...)
{
  .conductance_graph_impl(x = x, fit = fit, quantile = quantile, ...)
}

#for validation and debugging, generate a simple 1D lattice
Lattice1D <- function(spdat, coords, fn)
{
  stopifnot(all(!is.na(spdat)))
  coords_mat <- .coords_matrix(coords, NULL)
  stopifnot(all(coords_mat >= 0 & coords_mat <= 1))

  rl <- rast(nrows = 1, ncols = nrow(spdat), nlyrs = ncol(spdat),
             xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  values(rl) <- spdat
  names(rl) <- paste0("var", seq_len(ncol(spdat)))

  conductance_surface(rl, coords_mat)
}
