.validate_scale_vector <- function(x, names, arg)
{
  if (is.null(x))
    return(NULL)
  x <- as.numeric(x)
  if (length(x) == 1L)
    x <- rep(x, length(names))
  if (length(x) != length(names))
    stop("`", arg, "` must have length 1 or match the number of raster layers")
  names(x) <- names
  x
}

.validate_scale_grid <- function(scale_grid, names)
{
  if (is.null(scale_grid))
    return(NULL)
  if (!is.list(scale_grid))
    stop("`scale_grid` must be a named list of numeric candidate values")
  if (is.null(names(scale_grid)))
    names(scale_grid) <- names
  if (!setequal(names(scale_grid), names))
    stop("`scale_grid` names must match the raster layer names")
  scale_grid <- scale_grid[names]
  lapply(scale_grid, function(x)
  {
    x <- sort(unique(as.numeric(x)))
    if (!length(x) || anyNA(x))
      stop("Each `scale_grid` entry must contain at least one finite numeric value")
    x
  })
}

.resolve_scale_fun <- function(scale_fun, scale_args = NULL)
{
  scale_args <- if (is.null(scale_args)) list() else as.list(scale_args)
  if (!is.null(scale_fun))
    return(list(fun = scale_fun, args = scale_args, source = "user"))

  if (!requireNamespace("multiScaleR", quietly = TRUE))
    stop("Install package `multiScaleR` or supply a custom `scale_fun`.",
         call. = FALSE)

  kernel <- scale_args$kernel
  if (is.null(kernel))
    kernel <- "gaussian"
  scale_args$kernel <- NULL

  list(
    fun = function(r, scale, kernel = "gaussian", ...)
      multiScaleR::kernel_scale.raster(r, sigma = scale, kernel = kernel, ...),
    args = c(list(kernel = kernel), scale_args),
    source = "multiScaleR::kernel_scale.raster"
  )
}

.scale_covariate_stack <- function(covariates, scales, scale_fun, scale_args = NULL)
{
  covariates <- .as_spatraster(covariates)
  stopifnot(length(scales) == terra::nlyr(covariates))

  scaled <- lapply(seq_len(terra::nlyr(covariates)), function(i)
  {
    out <- do.call(scale_fun, c(list(covariates[[i]], scale = scales[[i]]), scale_args))
    out <- .as_spatraster(out)
    if (terra::nlyr(out) != 1L)
      stop("`scale_fun` must return a single-layer SpatRaster for each input layer")
    out
  })
  out <- do.call(c, scaled)
  names(out) <- names(covariates)
  out
}

.update_surface_covariates <- function(surface, covariates)
{
  covariates <- .as_spatraster(covariates)
  if (!all(names(covariates) == surface$covariates))
    stop("Scaled raster layer names must match `data$covariates`")

  xy <- surface$vertex_coordinates
  if (is.null(xy))
    stop("Surface does not retain `vertex_coordinates`, so covariates cannot be updated in place")

  cells <- terra::cellFromXY(covariates[[1]], xy)
  if (anyNA(cells))
    stop("Scaled raster no longer matches the original graph extent")

  spdat <- terra::values(covariates, dataframe = FALSE)
  x <- spdat[cells, , drop = FALSE]
  if (anyNA(x))
    stop("Scaled raster changed the missing-data pattern on active graph cells")

  is_factor <- vapply(seq_len(terra::nlyr(covariates)),
                      function(i) is.factor(covariates[[i]]),
                      logical(1))
  x <- as.data.frame(x)
  if (any(is_factor))
  {
    factors <- names(covariates)[is_factor]
    for (i in factors)
    {
      lev <- .factor_levels(covariates[[i]])
      x[, i] <- factor(lev$labels[match(x[, i], lev$ids)])
    }
  }

  out <- surface
  out$x <- x
  out$covariates <- colnames(x)
  if (!is.null(surface$stack))
    out$stack <- covariates
  out
}

.evaluate_scale_candidate <- function(scales,
                                      formula,
                                      response_name,
                                      response_value,
                                      covariates,
                                      coords,
                                      directions,
                                      conductance_model,
                                      measurement_model,
                                      objective,
                                      theta,
                                      postprocess,
                                      scale_fun,
                                      scale_args,
                                      terradish_args,
                                      cache_env,
                                      eval_env)
{
  key <- paste(formatC(scales, digits = 8, format = "fg", flag = "#"), collapse = "|")
  if (exists(key, envir = cache_env, inherits = FALSE))
    return(get(key, envir = cache_env, inherits = FALSE))

  scaled_covariates <- .scale_covariate_stack(covariates, scales, scale_fun, scale_args)
  if (!is.null(postprocess))
    scaled_covariates <- postprocess(scaled_covariates)

  surface <- tryCatch(
    .update_surface_covariates(eval_env$base_surface, scaled_covariates),
    error = function(e) NULL
  )
  if (is.null(surface))
    surface <- conductance_surface(scaled_covariates, coords,
                                   directions = directions, saveStack = TRUE)

  fit <- tryCatch(
    suppressWarnings(
      local({
        assign(response_name, response_value, envir = environment())
        do.call(
          terradish,
          c(list(
            formula = formula,
            data = surface,
            conductance_model = conductance_model,
            measurement_model = measurement_model,
            theta = theta
          ), terradish_args)
        )
      })
    ),
    error = function(e) e
  )

  result <- if (inherits(fit, "error"))
  {
    list(
      scales = scales,
      score = Inf,
      fit = NULL,
      surface = surface,
      status = paste("ERROR:", conditionMessage(fit)),
      aic = NA_real_,
      loglik = NA_real_
    )
  }
  else
  {
    score <- switch(objective,
                    aic = fit$aic,
                    logLik = -fit$loglik)
    list(
      scales = scales,
      score = score,
      fit = fit,
      surface = surface,
      status = "OK",
      aic = fit$aic,
      loglik = fit$loglik
    )
  }

  assign(key, result, envir = cache_env)
  result
}

.scale_evaluations_table <- function(results)
{
  if (!length(results))
    return(data.frame())
  scales <- do.call(rbind, lapply(results, `[[`, "scales"))
  out <- data.frame(
    do.call(cbind, lapply(seq_len(ncol(scales)), function(i) scales[, i])),
    objective = vapply(results, `[[`, numeric(1), "score"),
    aic = vapply(results, `[[`, numeric(1), "aic"),
    loglik = vapply(results, `[[`, numeric(1), "loglik"),
    status = vapply(results, `[[`, character(1), "status"),
    stringsAsFactors = FALSE
  )
  names(out)[seq_len(ncol(scales))] <- colnames(scales)
  out
}

#' Optimize raster scales of effect around a terradish fit
#'
#' Performs an outer optimization over raster scale parameters by repeatedly
#' rescaling raster covariates, rebuilding or updating the conductance surface,
#' and refitting \code{\link{terradish}}.
#'
#' @param formula A terradish model formula.
#' @param covariates A \code{terra::SpatRaster} containing the raster layers
#'   whose scales of effect will be optimized.
#' @param coords Focal coordinates passed to \code{\link{conductance_surface}}.
#' @param scales Optional starting values for the scale parameters. Required for
#'   \code{search = "coordinate"} unless \code{scale_grid} is supplied.
#' @param lower,upper Lower and upper bounds for each scale parameter. Used to
#'   construct \code{scale_grid} when not supplied, and required for
#'   \code{search = "coordinate"}.
#' @param scale_grid Optional named list of candidate scale values for
#'   \code{search = "grid"}.
#' @param search One of \code{"coordinate"} or \code{"grid"}.
#' @param objective Objective used for the outer optimization:
#'   \code{"aic"} or \code{"logLik"}.
#' @param grid_points Number of equally spaced candidate scales per raster when
#'   \code{search = "grid"} and \code{scale_grid} is not supplied.
#' @param maxit Maximum number of outer coordinate-search iterations.
#' @param tol Convergence tolerance for coordinate search, measured as the
#'   maximum absolute change in the scale vector between iterations.
#' @param scale_fun Optional function used to rescale one raster layer at a
#'   time. It must accept a single-layer raster as its first argument and a
#'   numeric \code{scale=} argument. If \code{NULL}, the function tries to use
#'   \code{multiScaleR::kernel_scale.raster()}.
#' @param scale_args Optional named list passed to \code{scale_fun}.
#' @param postprocess Optional function applied to the full scaled raster stack
#'   before constructing the conductance surface. Defaults to
#'   \code{\link{scale_covariates}} so each scaled layer is re-standardized.
#'   Set to \code{NULL} to use the scaled rasters as returned by
#'   \code{scale_fun}.
#' @param directions Neighborhood definition passed to
#'   \code{\link{conductance_surface}}.
#' @param theta Optional starting values passed to \code{\link{terradish}} at
#'   each outer evaluation.
#' @param verbose Print outer-optimization progress?
#' @param ... Additional arguments passed to \code{\link{terradish}}.
#'
#' @details
#' This helper treats scale optimization as an outer problem around the usual
#' \code{\link{terradish}} fit, rather than as an inner Newton/BFGS parameter.
#' That makes it compatible with raster-scaling approaches such as
#' \code{multiScaleR::kernel_scale.raster()}, while keeping the terradish core
#' optimizer unchanged.
#'
#' Coordinate search updates one raster scale at a time using
#' \code{\link[stats]{optimize}} while holding the others fixed. Grid search
#' evaluates every combination in \code{scale_grid}. In both modes, failed
#' terradish fits receive an infinite outer objective and are retained in the
#' returned evaluation table.
#'
#' @return A list with components including:
#' \item{par}{The best scale vector found.}
#' \item{value}{The best outer objective value.}
#' \item{fit}{The best terradish fit.}
#' \item{surface}{The conductance surface used by the best fit.}
#' \item{evaluations}{A data frame of all outer evaluations.}
#' \item{search}{The search strategy used.}
#' \item{objective}{The outer objective used for ranking candidate scales.}
#'
#' @examples
#' \dontrun{
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#' covariates <- c(melip.altitude, melip.forestcover)
#' names(covariates) <- c("altitude", "forestcover")
#'
#' fit_scale <- terradish_scale_optim(
#'   melip.Fst ~ altitude + forestcover,
#'   covariates = covariates,
#'   coords = melip.coords,
#'   lower = c(altitude = 50, forestcover = 50),
#'   upper = c(altitude = 500, forestcover = 500),
#'   search = "grid",
#'   grid_points = 3,
#'   measurement_model = mlpe,
#'   optimizer = "bfgs"
#' )
#' fit_scale$par
#' }
#'
#' @export
terradish_scale_optim <- function(formula,
                                  covariates,
                                  coords,
                                  scales = NULL,
                                  lower = NULL,
                                  upper = NULL,
                                  scale_grid = NULL,
                                  search = c("coordinate", "grid"),
                                  objective = c("aic", "logLik"),
                                  grid_points = 5L,
                                  maxit = 10L,
                                  tol = 0.1,
                                  scale_fun = NULL,
                                  scale_args = NULL,
                                  postprocess = scale_covariates,
                                  directions = 8L,
                                  theta = NULL,
                                  verbose = TRUE,
                                  ...)
{
  covariates <- .as_spatraster(covariates)
  search <- match.arg(search)
  objective <- match.arg(objective)
  terms_obj <- terms(formula)
  response_idx <- attr(terms_obj, "response")
  if (!response_idx)
    stop("`formula` must have a response matrix on the left-hand side")
  response_name <- all.vars(formula[[2]])
  if (length(response_name) != 1L)
    stop("The left-hand side of `formula` must be a single matrix object")
  response_env <- environment(formula)
  if (is.null(response_env))
    response_env <- parent.frame()
  response_value <- get(response_name, envir = response_env, inherits = TRUE)
  scale_names <- names(covariates)
  if (is.null(scale_names))
    stop("`covariates` must have named raster layers")

  scales <- .validate_scale_vector(scales, scale_names, "scales")
  lower <- .validate_scale_vector(lower, scale_names, "lower")
  upper <- .validate_scale_vector(upper, scale_names, "upper")
  scale_grid <- .validate_scale_grid(scale_grid, scale_names)
  grid_points <- as.integer(grid_points)[1]
  maxit <- as.integer(maxit)[1]
  stopifnot(!is.na(grid_points), grid_points >= 2L)
  stopifnot(!is.na(maxit), maxit >= 1L)
  stopifnot(is.numeric(tol), length(tol) == 1L, is.finite(tol), tol >= 0)

  if (identical(search, "coordinate"))
  {
    if (is.null(scales))
    {
      if (is.null(lower) || is.null(upper))
        stop("`scales`, or both `lower` and `upper`, are required for coordinate search")
      scales <- (lower + upper) / 2
    }
    if (is.null(lower) || is.null(upper))
      stop("`lower` and `upper` are required for coordinate search")
  }

  if (identical(search, "grid") && is.null(scale_grid))
  {
    if (is.null(lower) || is.null(upper))
      stop("Supply `scale_grid`, or both `lower` and `upper`, for grid search")
    scale_grid <- stats::setNames(
      lapply(seq_along(scale_names), function(i)
        seq(lower[[i]], upper[[i]], length.out = grid_points)),
      scale_names
    )
  }

  resolved_scale_fun <- .resolve_scale_fun(scale_fun, scale_args)
  terradish_args <- list(...)
  conductance_model <- if (is.null(terradish_args$conductance_model))
    loglinear_conductance
  else
    terradish_args$conductance_model
  measurement_model <- if (is.null(terradish_args$measurement_model))
    mlpe
  else
    terradish_args$measurement_model
  terradish_args$conductance_model <- NULL
  terradish_args$measurement_model <- NULL
  terradish_args$theta <- NULL
  base_covariates <- if (!is.null(postprocess)) postprocess(covariates) else covariates
  base_surface <- conductance_surface(base_covariates, coords,
                                      directions = directions, saveStack = TRUE)
  cache_env <- new.env(parent = emptyenv())
  eval_env <- new.env(parent = emptyenv())
  eval_env$base_surface <- base_surface
  all_results <- list()

  evaluate_candidate <- function(scale_values)
  {
    scale_values <- .validate_scale_vector(scale_values, scale_names, "scale_values")
    result <- .evaluate_scale_candidate(
      scales = scale_values,
      formula = formula,
      response_name = response_name,
      response_value = response_value,
      covariates = covariates,
      coords = coords,
      directions = directions,
      conductance_model = conductance_model,
      measurement_model = measurement_model,
      objective = objective,
      theta = theta,
      postprocess = postprocess,
      scale_fun = resolved_scale_fun$fun,
      scale_args = resolved_scale_fun$args,
      terradish_args = terradish_args,
      cache_env = cache_env,
      eval_env = eval_env
    )
    all_results[[length(all_results) + 1L]] <<- result
    if (isTRUE(verbose))
      cat("scales:", paste(sprintf("%s=%.3f", names(scale_values), scale_values), collapse = ", "),
          "| objective:", sprintf("%.4f", result$score),
          "| status:", result$status, "\n")
    result
  }

  best <- NULL
  iter <- 0L
  converged <- FALSE

  if (identical(search, "grid"))
  {
    combos <- expand.grid(scale_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    combos <- as.matrix(combos[, scale_names, drop = FALSE])
    for (i in seq_len(nrow(combos)))
    {
      result <- evaluate_candidate(combos[i, ])
      if (is.null(best) || result$score < best$score)
        best <- result
    }
    iter <- nrow(combos)
    converged <- TRUE
  }
  else
  {
    current <- scales
    best <- evaluate_candidate(current)
    for (iter_i in seq_len(maxit))
    {
      previous <- current
      for (j in seq_along(current))
      {
        fn_j <- function(x)
          evaluate_candidate(replace(current, j, x))$score
        opt_j <- stats::optimize(fn_j,
                                 interval = c(lower[[j]], upper[[j]]))
        current[[j]] <- opt_j$minimum
      }
      candidate <- evaluate_candidate(current)
      if (candidate$score < best$score)
        best <- candidate
      iter <- iter_i
      if (max(abs(current - previous)) <= tol)
      {
        converged <- TRUE
        break
      }
    }
  }

  evaluations <- unique(.scale_evaluations_table(all_results))
  out <- list(
    call = match.call(),
    par = best$scales,
    value = best$score,
    fit = best$fit,
    surface = best$surface,
    evaluations = evaluations,
    search = search,
    objective = objective,
    scale_fun = resolved_scale_fun$source,
    iterations = iter,
    convergence = if (isTRUE(converged)) 0L else 1L
  )
  class(out) <- "terradish_scale_optim"
  out
}

#' @rdname terradish_scale_optim
#' @export
terra_radish_scale_optim <- terradish_scale_optim
