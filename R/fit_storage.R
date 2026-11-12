#' Shrink a fitted model for storage
#'
#' Removes components that are useful for prediction or diagnostics but can
#' make a fitted \code{terradish} object expensive to serialize. The fitted
#' coefficients, likelihood, Hessian, fitted response, convergence diagnostics,
#' and approximation metadata are retained unless explicitly dropped.
#'
#' @param object A fitted \code{terradish} object.
#' @param drop Character vector naming optional components to remove. Supported
#'   values are \code{"submodels"}, \code{"response"}, and
#'   \code{"leverage"}. Dropping \code{"submodels"} removes the fitted
#'   conductance and measurement-model closures, which usually account for most
#'   of the serialized object size. Dropping \code{"response"} also removes the
#'   response copy stored in the comparison contract.
#'
#' @details
#' A slim object remains suitable for \code{print()}, \code{summary()},
#' \code{coef()}, \code{vcov()}, \code{confint()}, \code{logLik()}, \code{AIC()}, and fitted-value
#' extraction and residual-permutation simulation. Operations that must
#' reevaluate the conductance or measurement model, including conductance-
#' surface prediction and most plot types, require the retained model closures
#' and will stop with an informative error if those closures were dropped.
#'
#' Slimming is irreversible. Retain the original fit, model formula, graph, and
#' response when later prediction, simulation, or refitting may be needed.
#'
#' @return A fitted \code{terradish} object with a \code{storage} element that
#'   records the removed components and the object's size before and after
#'   slimming.
#'
#' @seealso \code{\link{terradish}} (argument \code{slim}),
#'   \code{\link{terradish_inference}}
#'
#' @family solver settings
#'
#' @examples
#' \donttest{
#' data(melip)
#' covariates <- c(terra::unwrap(melip.altitude),
#'                 terra::unwrap(melip.forestcover))
#' coords <- terra::unwrap(melip.coords)
#' names(covariates) <- c("altitude", "forestcover")
#' # Coarsened so the example runs quickly; name layers, then scale.
#' covariates <- scale_covariates(terra::aggregate(covariates, fact = 3,
#'                                                 na.rm = TRUE))
#' # Keep one site per raster cell (some melip sites share a cell).
#' keep <- !duplicated(terra::cellFromXY(covariates, terra::crds(coords)))
#' coords <- coords[keep]
#' Fst <- melip.Fst[keep, keep]
#' surface <- conductance_surface(covariates, coords)
#' fit <- terradish(Fst ~ altitude + forestcover, surface,
#'                  measurement_model = mlpe)
#' fit_slim <- slim_terradish(fit)
#' fit_slim$storage[c("bytes_before", "bytes_after")]
#' coef(fit_slim)     # estimates, summaries, and intervals still work
#' confint(fit_slim)
#' }
#'
#' @export
slim_terradish <- function(object,
                           drop = c("submodels", "leverage"))
{
  if (!inherits(object, c("terradish", "radish")))
    stop("`object` must be a fitted terradish model.", call. = FALSE)

  allowed <- c("submodels", "response", "leverage")
  drop <- unique(as.character(drop))
  unknown <- setdiff(drop, allowed)
  if (length(unknown))
    stop("Unknown component in `drop`: ", paste(unknown, collapse = ", "),
         call. = FALSE)

  before <- as.numeric(utils::object.size(object))
  if ("submodels" %in% drop)
    object$submodels <- NULL
  if ("leverage" %in% drop)
    object$leverage <- NULL
  if ("response" %in% drop)
  {
    object$response <- NULL
    object$fit$response <- NULL
    object$comparison$response <- NULL
  }

  object$storage <- list(
    slim = TRUE,
    dropped = drop,
    bytes_before = before,
    bytes_after = NA_real_
  )
  object$storage$bytes_after <- as.numeric(utils::object.size(object))
  object
}

.terradish_has_submodels <- function(object)
{
  is.list(object$submodels) &&
    is.function(object$submodels$f) &&
    is.function(object$submodels$g)
}

.terradish_require_submodels <- function(object, operation)
{
  if (!.terradish_has_submodels(object))
    stop(operation, " requires model closures that were removed by `slim_terradish()`. ",
         "Use the unslimmed fit for this operation.", call. = FALSE)
  invisible(TRUE)
}
