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
#' Reads spatial covariates and focal-point coordinates and builds a weighted
#' graph suitable for computing conductance-weighted resistance distances.
#' This is the first step of every \code{terradish} workflow.
#'
#' @param covariates A \code{SpatRaster} (or \code{PackedSpatRaster}) whose
#'   layers contain the spatial covariates to be used in the conductance model.
#'   Layer names must match the right-hand side of the formula passed to
#'   \code{\link{terradish}}.  Scale continuous layers with
#'   \code{\link{scale_covariates}} before calling this function.
#' @param coords Locations of the focal sampling sites, in the same coordinate
#'   reference system as \code{covariates}.  Supported inputs:
#'   \code{terra::SpatVector} of points, a two-column numeric matrix (x in
#'   column 1, y in column 2), or a data frame with columns named \code{x} and
#'   \code{y} (or at least two columns, using the first two).
#' @param directions Adjacency rule for the graph.  \code{4} connects each
#'   cell to its four horizontal/vertical neighbours (rook adjacency);
#'   \code{8} additionally includes the four diagonal neighbours (queen
#'   adjacency).  \strong{\code{8} is recommended} for landscape genetics:
#'   it allows diagonal movement, produces smoother resistance surfaces, and
#'   is less sensitive to grid orientation artifacts.
#' @param saveStack Logical.  If \code{TRUE} (default), the masked
#'   \code{SpatRaster} is stored inside the returned \code{terradish_graph}
#'   object.  Required for raster-based conductance visualization
#'   (\code{plot(fit, type = "surface")}) and marginal-effect plots.  Set to
#'   \code{FALSE} to reduce memory use when the raster is large and
#'   visualization is not needed.
#' @param crop_buffer Optional nonnegative map-unit buffer to crop
#'   \code{covariates} around the focal coordinates before graph construction.
#'   When supplied, \code{\link{crop_to_focal_buffer}} is called internally.
#'   A scalar applies the same buffer in both x and y directions; a
#'   length-two vector gives \code{c(x_buffer, y_buffer)}.  Units are the
#'   raster's map units (degrees for lon/lat, metres for projected rasters).
#'   See Details.
#'
#' @details
#' \strong{Missing values:} \code{NA} is propagated uniformly across all
#' layers so every cell is either complete or entirely missing.  A warning is
#' issued if some cells are \code{NA} in only a subset of layers, because this
#' indicates that different covariates cover different areas — model comparison
#' across raster stacks with different \code{NA} patterns is problematic since
#' the graph vertices differ.  Always compare models using the same
#' \code{terradish_graph} object.
#'
#' \strong{Disconnected components:} after masking, the graph may contain
#' cells that are isolated from the main component.  The function retains only
#' the largest connected component and warns about the number of pruned cells.
#' The function aborts if any focal site belongs to a disconnected component.
#'
#' \strong{Cropping:} when \code{crop_buffer} is supplied, cropping changes
#' the graph domain.  Always use the same \code{terradish_graph} for all
#' models you intend to compare.  A buffer that is too tight may cut off
#' landscape corridors that affect resistance distances; check stability by
#' refitting with several buffer values.
#'
#' \strong{Categorical covariates:} raster layers representing categorical
#' variables must be stored as factor-valued \code{SpatRaster} layers (see
#' \code{\link[terra]{as.factor}} and the example below).  Level names are
#' taken from the \code{VALUE} column of the levels table when present;
#' otherwise the first non-\code{ID} column is used.  Factor levels are
#' dummy-coded by \code{\link[stats]{model.matrix}} during model fitting, with
#' one level as the reference.
#'
#' \strong{Directions:} \code{directions = 4} gives a more conservative graph
#' that can underestimate resistance in diagonal corridors.
#' \code{directions = 8} is preferred unless there is a specific reason to
#' restrict movement to the cardinal directions.
#'
#' @return An object of class \code{"terradish_graph"} (also inheriting
#'   \code{"radish_graph"}) containing:
#' \describe{
#'   \item{\code{demes}}{Integer indices of cells corresponding to focal sites.}
#'   \item{\code{x}}{Data frame of covariate values at non-missing cells.}
#'   \item{\code{adj}}{Edge list (upper-triangular, 0-based).}
#'   \item{\code{edge_pairs}}{Integer matrix of adjacent cell pairs.}
#'   \item{\code{vertex_coordinates}}{xy-coordinate matrix of active cells.}
#'   \item{\code{directions}}{The adjacency rule (4 or 8).}
#'   \item{\code{covariates}}{Names of the covariate columns.}
#'   \item{\code{stack}}{The masked \code{SpatRaster} (or \code{NULL} if
#'     \code{saveStack = FALSE}).}
#' }
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
#' Extracts per-cell conductance values and pointwise confidence intervals from
#' a fitted \code{terradish} model, returning either a raster (when the graph
#' stores its stack) or a matrix.
#'
#' @param x A \code{terradish_graph} object returned by
#'   \code{\link{conductance_surface}}.
#' @param fit A fitted model returned by \code{\link{terradish}}.
#' @param quantile Confidence level for the pointwise conductance interval.
#'   The default \code{0.95} gives a 95% interval computed via the delta method
#'   applied to the MLE of \code{theta} and its inverse Hessian.
#' @param support Optional support constraint for conductance evaluation.
#'   Use \code{"none"} (default) to evaluate the fitted conductance model over
#'   the full graph covariate domain, or \code{"focal"} to clamp selected
#'   covariates to focal-site support.
#' @param support_probs Optional length-2 probability vector used with
#'   \code{support = "focal"} to compute clamping bounds from focal-cell
#'   quantiles. The default \code{c(0, 1)} uses the focal min/max.
#' @param clamp_covariates Optional character vector of covariate names to
#'   clamp when \code{support = "focal"}. Defaults to all numeric covariates in
#'   \code{x$x}.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Confidence intervals are computed on the log-conductance scale and then
#' exponentiated, so the point estimate always lies within the interval.  The
#' intervals are asymptotic (delta method) and reflect only uncertainty in the
#' conductance parameters \code{theta}; they do not account for uncertainty in
#' the nuisance parameters \code{phi}.
#'
#' With \code{support = "focal"}, selected covariates are clamped to focal-site
#' support before conductance is evaluated. This constrains prediction to
#' empirical covariate support among \code{x$demes} and can stabilize extreme
#' tails when focal sampling is sparse.
#'
#' The function requires the MLE to be in the interior of the parameter space
#' (\code{fit$fit$boundary == FALSE}); if the MLE is on the boundary (no
#' detectable IBR signal) an error is thrown.
#'
#' @return
#' If \code{x$stack} is non-\code{NULL} (the default when
#' \code{saveStack = TRUE} in \code{\link{conductance_surface}}), a three-layer
#' \code{terra::SpatRaster} with layer names:
#' \describe{
#'   \item{\code{est}}{Point-estimate conductance at each active cell.}
#'   \item{\code{lower<Q>}}{Lower bound of the \code{quantile*100}\% interval
#'     (e.g. \code{lower95} for the default \code{quantile = 0.95}).}
#'   \item{\code{upper<Q>}}{Upper bound of the interval.}
#' }
#' Otherwise, a numeric matrix with one row per active cell and columns
#' \code{est}, \code{lower<Q>}, and \code{upper<Q>}.
#'
#' @seealso \code{\link{terradish}}, \code{\link{conductance_surface}},
#'   \code{\link{plot.terradish}}
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
#' names(cond)   # "est", "lower95", "upper95"
#' plot(cond)
#'
#' @export
conductance <- function(x, ...)
{
  UseMethod("conductance")
}

.conductance_graph_impl <- function(x, fit, quantile = 0.95,
                                    support = c("none", "focal"),
                                    support_probs = c(0, 1),
                                    clamp_covariates = NULL,
                                    ...)
{
  support <- match.arg(support)
  stopifnot(inherits(fit, c("terradish", "radish")))
  stopifnot(!fit$fit$boundary && !is.null(fit$mle$theta))

  graph_eval <- .clamp_graph_covariates(
    data = x,
    support = support,
    support_probs = support_probs,
    clamp_covariates = clamp_covariates
  )

  conductance_model <- fit$submodels$f
  if (!identical(support, "none"))
  {
    conductance_factory <- fit$submodels$f_factory
    if (!inherits(conductance_factory,
                  c("terradish_conductance_model_factory",
                    "radish_conductance_model_factory")))
      stop("This `fit` does not store a reusable conductance-model factory for support-constrained prediction.",
           " Refit the model with the current terradish version and retry.",
           call. = FALSE)

    rebuilt_internal <- conductance_factory(fit$formula, graph_eval$x)
    conductance_model <- .externalize_conductance_model(rebuilt_internal)
  }

  conductance <- conductance_model(fit$mle$theta)
  ci <- conductance$confint(theta = fit$mle$theta, vcov = -solve(fit$mle$hessian), 
                            quantile = quantile, scale = "conductance")
  colnames(ci) <- paste0(c("lower", "upper"), round(100*quantile, 1))

  if (!is.null(graph_eval$stack))
  {
    template <- graph_eval$stack[[1]]
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
conductance.terradish_graph <- function(x, fit, quantile = 0.95,
                                        support = c("none", "focal"),
                                        support_probs = c(0, 1),
                                        clamp_covariates = NULL,
                                        ...)
{
  .conductance_graph_impl(x = x, fit = fit, quantile = quantile,
                          support = support,
                          support_probs = support_probs,
                          clamp_covariates = clamp_covariates,
                          ...)
}

#' @rdname conductance
#' @method conductance radish_graph
#' @export
conductance.radish_graph <- function(x, fit, quantile = 0.95,
                                     support = c("none", "focal"),
                                     support_probs = c(0, 1),
                                     clamp_covariates = NULL,
                                     ...)
{
  .conductance_graph_impl(x = x, fit = fit, quantile = quantile,
                          support = support,
                          support_probs = support_probs,
                          clamp_covariates = clamp_covariates,
                          ...)
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
