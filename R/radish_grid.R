# Internal worker helpers shared by serial and parallel paths.
.terradish_grid_worker <- function(idx, state)
{
  worker_formula <- reformulate(state$term_labels)
  conductance_model_local <- if (!is.null(state$conductance_model_factory))
    state$conductance_model_factory(worker_formula, state$data$x)
  else
    state$conductance_model
  phi_start <- NULL
  out <- vector("list", length(idx))
  for (j in seq_along(idx))
  {
    i <- idx[j]
    out[[j]] <- tryCatch({
      alg_args <- list(f = conductance_model_local,
                       g = state$measurement_model,
                       s = state$data,
                       S = state$S,
                       theta = c(state$theta[i,]),
                       nu = state$nu,
                       gradient = FALSE,
                       hessian = FALSE,
                       partial = FALSE,
                       nonnegative = state$nonnegative)
      if (!is.null(phi_start))
        alg_args$phi <- phi_start
      obj <- do.call(state$radish_algorithm, alg_args)
      phi_start <- obj$phi
      list(i = i,
           loglik = -obj$objective,
           phi = obj$phi,
           covariance = if (state$covariance) as.matrix(obj$covariance) else NULL)
    }, error = function(e) list(i = i, error = paste(conditionMessage(e), " @ ", deparse(conditionCall(e)))))
  }
  out
}

.terradish_distance_worker <- function(idx, state)
{
  worker_formula <- reformulate(state$term_labels)
  use_namespace_workers <- isTRUE(state$use_namespace_workers)
  conductance_model_factory <- state$conductance_model_factory
  if (!is.null(state$conductance_model_factory_name))
  {
    if (use_namespace_workers)
      conductance_model_factory <- get(state$conductance_model_factory_name,
                                       envir = asNamespace("terradish"),
                                       inherits = FALSE)
    else
      conductance_model_factory <- get(state$conductance_model_factory_name,
                                       envir = globalenv(),
                                       mode = "function",
                                       inherits = TRUE)
  }
  conductance_model_local <- if (!is.null(conductance_model_factory))
    conductance_model_factory(worker_formula, state$data$x)
  else
    state$conductance_model
  out <- vector("list", length(idx))
  for (j in seq_along(idx))
  {
    i <- idx[j]
    out[[j]] <- tryCatch({
      obj <- state$radish_algorithm(f = conductance_model_local,
                                    g = state$leastsquares,
                                    s = state$data,
                                    S = diag(length(state$data$demes)),
                                    theta = c(state$theta[i,]),
                                    objective = FALSE,
                                    gradient = FALSE,
                                    hessian = FALSE,
                                    partial = FALSE,
                                    nonnegative = TRUE)
      cov_i <- as.matrix(obj$covariance)
      list(i = i, value = if (state$covariance) cov_i else state$dist_from_cov(cov_i))
    }, error = function(e) list(i = i, error = paste(conditionMessage(e), " @ ", deparse(conditionCall(e)))))
  }
  out
}

.parallel_fun_name <- function(fun, candidates, namespace = "terradish")
{
  ns_env <- asNamespace(namespace)
  for (nm in candidates)
  {
    obj <- get0(nm, envir = ns_env, mode = "function", inherits = FALSE)
    if (!is.null(obj) &&
        identical(formals(fun), formals(obj)) &&
        identical(body(fun), body(obj)))
      return(nm)
  }
  NULL
}

.resolve_parallel_fun <- function(fun, fun_name, use_namespace_workers, namespace = "terradish")
{
  if (!is.null(fun_name))
  {
    if (use_namespace_workers)
      return(get(fun_name, envir = asNamespace(namespace), inherits = FALSE))
    return(get(fun_name, mode = "function", inherits = TRUE))
  }
  fun
}

.use_namespace_workers <- function(namespace = "terradish")
{
  if (file.exists(file.path(getwd(), "R", "radish_grid.R")))
    return(FALSE)

  if (!(namespace %in% loadedNamespaces()))
    return(FALSE)

  ns_path <- tryCatch(getNamespaceInfo(asNamespace(namespace), "path"),
                      error = function(e) NULL)
  if (is.null(ns_path))
    return(FALSE)

  cwd <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE),
                  error = function(e) NULL)
  ns_path <- tryCatch(normalizePath(ns_path, winslash = "/", mustWork = FALSE),
                      error = function(e) ns_path)

  !identical(ns_path, cwd)
}

.validate_theta_grid <- function(theta, parameter_names)
{
  stopifnot(is.matrix(theta))
  stopifnot(length(parameter_names) == ncol(theta))

  theta_names <- colnames(theta)
  if (is.null(theta_names))
  {
    colnames(theta) <- parameter_names
    return(theta)
  }

  if (anyDuplicated(theta_names))
    stop("`theta` column names must be unique")

  missing_names <- setdiff(parameter_names, theta_names)
  extra_names <- setdiff(theta_names, parameter_names)
  if (length(missing_names) || length(extra_names))
    stop("`theta` column names must match the conductance-model parameters")

  theta[, parameter_names, drop = FALSE]
}

.terradish_grid_approximation <- function(data, S, approximation = c("none", "landmark", "coarse_raster"),
                                          approximation_control = NULL, covariance = FALSE)
{
  approximation <- match.arg(approximation)
  if (identical(approximation, "none"))
    return(list(data = data,
                S = S,
                info = list(type = "none", used = FALSE)))

  if (identical(approximation, "landmark"))
  {
    subset <- .terradish_landmark_subset(data, S, approximation_control = approximation_control)
    info <- list(
      type = "landmark",
      used = isTRUE(subset$used),
      n_landmarks = subset$control$n_landmarks,
      full_focal = nrow(S),
      focal_fraction = subset$control$n_landmarks / nrow(S),
      method = subset$control$method,
      landmark_index = subset$index
    )

    if (isTRUE(covariance) && isTRUE(subset$used))
      warning("`covariance = TRUE` with `approximation = \"landmark\"` returns landmark covariance only.",
              call. = FALSE)

    return(list(data = subset$data,
                S = subset$S,
                info = info))
  }

  coarse <- .coarse_raster_surface(data, approximation_control = approximation_control)
  info <- list(
    type = "coarse_raster",
    used = isTRUE(coarse$used),
    factor = coarse$control$factor,
    directions = coarse$control$directions,
    full_vertices = coarse$full_vertices,
    coarse_vertices = coarse$coarse_vertices,
    vertex_fraction = coarse$coarse_vertices / coarse$full_vertices,
    full_focal = nrow(S),
    unique_demes = coarse$unique_demes,
    duplicate_demes = coarse$duplicate_demes
  )

  if (isTRUE(covariance) && isTRUE(coarse$used))
    warning("`covariance = TRUE` with `approximation = \"coarse_raster\"` returns covariance on the aggregated graph.",
            call. = FALSE)

  list(data = coarse$surface,
       S = S,
       info = info)
}

.terradish_grid_parallel_chunk <- function(idx, theta, term_labels, data,
                                           conductance_model_factory_name = NULL,
                                           conductance_model_factory = NULL,
                                           measurement_model_name = NULL,
                                           measurement_model = NULL,
                                           S, nu, nonnegative, covariance,
                                           use_namespace_workers)
{
  resolve_fun <- function(fun, fun_name)
  {
    if (!is.null(fun_name))
    {
      if (use_namespace_workers)
        return(get(fun_name, envir = asNamespace("terradish"), inherits = FALSE))
      return(get(fun_name, envir = globalenv(), mode = "function", inherits = TRUE))
    }
    fun
  }
  formula <- reformulate(term_labels)
  conductance_model_factory <- resolve_fun(conductance_model_factory,
                                           conductance_model_factory_name)
  measurement_model <- resolve_fun(measurement_model,
                                   measurement_model_name)
  radish_algorithm <- if (use_namespace_workers)
    get("terradish_algorithm", envir = asNamespace("terradish"), inherits = FALSE)
  else
    get("terradish_algorithm", mode = "function", inherits = TRUE)
  conductance_model_local <- conductance_model_factory(formula, data$x)
  phi_start <- NULL
  out <- vector("list", length(idx))
  for (j in seq_along(idx))
  {
    i <- idx[j]
    out[[j]] <- tryCatch({
      alg_args <- list(f = conductance_model_local,
                       g = measurement_model,
                       s = data,
                       S = S,
                       theta = c(theta[i,]),
                       nu = nu,
                       gradient = FALSE,
                       hessian = FALSE,
                       partial = FALSE,
                       nonnegative = nonnegative)
      if (!is.null(phi_start))
        alg_args$phi <- phi_start
      obj <- do.call(terradish_algorithm, alg_args)
      phi_start <- obj$phi
      list(i = i,
           loglik = -obj$objective,
           phi = obj$phi,
           covariance = if (covariance) as.matrix(obj$covariance) else NULL)
    }, error = function(e) list(i = i, error = paste(conditionMessage(e), " @ ", deparse(conditionCall(e)))))
  }
  out
}

.terradish_grid_parallel_results <- function(theta, term_labels, data,
                                             conductance_model_factory, measurement_model,
                                             S, nu, nonnegative, covariance,
                                             cores,
                                             use_namespace_workers,
                                             worker_libpaths = .libPaths())
{
  idx <- seq_len(nrow(theta))
  splits <- split(idx, cut(idx, breaks = min(as.integer(cores), length(idx)), labels = FALSE))
  cl <- makeCluster(length(splits))
  on.exit(stopCluster(cl), add = TRUE)

  conductance_model_factory_name <- .parallel_fun_name(
    conductance_model_factory,
    c("loglinear_conductance", "linear_conductance")
  )
  measurement_model_name <- .parallel_fun_name(
    measurement_model,
    c("leastsquares", "mlpe", "generalized_wishart", "wishart_covariance")
  )

  if (use_namespace_workers)
  {
    clusterExport(cl, varlist = c("worker_libpaths"), envir = environment())
    clusterEvalQ(cl, {
      .libPaths(worker_libpaths)
      library(terradish)
      NULL
    })
  }
  else
  {
    worker_files <- c("R/backtracking.R", "R/hager_zhang.R", "R/newton_raphson.R",
                      "R/radish_conductance_model.R", "R/leastsquares.R",
                      "R/mlpe.R", "R/generalized_wishart.R", "R/wishart_covariance.R",
                      "R/radish_subproblem.R", "R/radish_algorithm.R")
    worker_wd <- getwd()
    clusterExport(cl, varlist = c("worker_files", "worker_wd"), envir = environment())
    clusterEvalQ(cl, {
      setwd(worker_wd)
      for (f in worker_files)
        source(f)
      NULL
    })
  }

  unlist(parLapply(
    cl, splits, .terradish_grid_parallel_chunk,
    theta = theta,
    term_labels = term_labels,
    data = data,
    conductance_model_factory_name = conductance_model_factory_name,
    conductance_model_factory = if (is.null(conductance_model_factory_name)) conductance_model_factory else NULL,
    measurement_model_name = measurement_model_name,
    measurement_model = if (is.null(measurement_model_name)) measurement_model else NULL,
    S = S,
    nu = nu,
    nonnegative = nonnegative,
    covariance = covariance,
    use_namespace_workers = use_namespace_workers
  ), recursive = FALSE)
}

.terradish_grid_serial_results <- function(idx, theta, model_formula, data,
                                           conductance_model_factory, measurement_model,
                                           S, nu, nonnegative, covariance)
{
  conductance_model_local <- conductance_model_factory(model_formula, data$x)
  phi_start <- NULL
  out <- vector("list", length(idx))
  for (j in seq_along(idx))
  {
    i <- idx[j]
    out[[j]] <- tryCatch({
      alg_args <- list(f = conductance_model_local,
                       g = measurement_model,
                       s = data,
                       S = S,
                       theta = c(theta[i,]),
                       nu = nu,
                       gradient = FALSE,
                       hessian = FALSE,
                       partial = FALSE,
                       nonnegative = nonnegative)
      if (!is.null(phi_start))
        alg_args$phi <- phi_start
      obj <- do.call(terradish_algorithm, alg_args)
      phi_start <- obj$phi
      list(i = i,
           loglik = -obj$objective,
           phi = obj$phi,
           covariance = if (covariance) as.matrix(obj$covariance) else NULL)
    }, error = function(e) list(i = i, error = paste(conditionMessage(e), " @ ", deparse(conditionCall(e)))))
  }
  out
}

.terradish_grid_impl <- function(theta, term_labels, data, S,
                                 conductance_model_factory, measurement_model,
                                 nu, nonnegative, covariance, cores,
                                 approximation = c("none", "landmark", "coarse_raster"),
                                 approximation_control = NULL)
{
  model_formula <- reformulate(term_labels)
  conductance_model <- conductance_model_factory(model_formula, data$x)
  default <- attr(conductance_model, "default")
  rm(conductance_model)

  stopifnot(ncol(theta) == length(default))
  theta <- .validate_theta_grid(theta, names(default))

  approx <- .terradish_grid_approximation(
    data = data,
    S = S,
    approximation = approximation,
    approximation_control = approximation_control,
    covariance = covariance
  )
  data_eval <- approx$data
  S_eval <- approx$S

  phi_template <- measurement_model(S = S_eval, E = rWishart(1, nrow(S_eval), diag(nrow(S_eval)))[,,1])$phi
  ll <- rep(NA_real_, nrow(theta))
  phi <- matrix(NA_real_, length(phi_template), nrow(theta))
  if (covariance)
    cv <- array(NA_real_, c(length(data_eval$demes), length(data_eval$demes), nrow(theta)))

  idx <- seq_len(nrow(theta))
  if (as.integer(cores) <= 1L || nrow(theta) <= 1L)
  {
    results <- .terradish_grid_serial_results(
      idx = idx,
      theta = theta,
      model_formula = model_formula,
      data = data_eval,
      conductance_model_factory = conductance_model_factory,
      measurement_model = measurement_model,
      S = S_eval,
      nu = nu,
      nonnegative = nonnegative,
      covariance = covariance
    )
  }
  else
  {
    use_namespace_workers <- .use_namespace_workers()
    results <- .terradish_grid_parallel_results(
      theta = theta,
      term_labels = term_labels,
      data = data_eval,
      conductance_model_factory = conductance_model_factory,
      measurement_model = measurement_model,
      S = S_eval,
      nu = nu,
      nonnegative = nonnegative,
      covariance = covariance,
      cores = cores,
      use_namespace_workers = use_namespace_workers,
      worker_libpaths = .libPaths()
    )
  }

  for (res in results)
  {
    if (!is.null(res$error))
      next
    ll[res$i] <- res$loglik
    phi[, res$i] <- res$phi
    if (covariance)
      cv[, , res$i] <- res$covariance
  }
  worker_errors <- unique(unlist(lapply(results, function(res) {
    if (is.null(res$error)) NA_character_ else res$error
  })))
  worker_errors <- worker_errors[!is.na(worker_errors)]
  if (length(worker_errors))
    warning("Parallel worker errors: ", paste(worker_errors, collapse = "; "))

  df <- list(theta = data.frame(theta),
             loglik = ll,
             phi = phi,
             covariance = if(!covariance) NULL else cv,
             approximation = approx$info)
  class(df) <- c("terradish_grid", "radish_grid")
  df
}

#' Evaluate likelihood of a parameterized conductance surface
#'
#' Calculates the profile likelihood of a parameterized conductance surface across
#' a grid of parameter values (e.g. the nuisance parameters are optimized at each
#' point on the grid).
#'
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param formula A formula with the name of a matrix of observed genetic distances on the lhs, and covariates in the creation of \code{data} on the rhs
#' @param data An object of class \code{terradish_graph} (see
#'   \code{\link{conductance_surface}})
#' @param conductance_model A function of class
#'   \code{terradish_conductance_model_factory} (see
#'   \code{\link{terradish_conductance_model_factory}})
#' @param measurement_model A function of class
#'   \code{terradish_measurement_model} (see
#'   \code{\link{terradish_measurement_model}})
#' @param nu Number of genetic markers (potentially used by \code{measurement_model})
#' @param nonnegative Force regression-like \code{measurement_model} to have nonnegative slope?
#' @param conductance Retained for backward compatibility. Only
#'   \code{conductance = TRUE} is currently implemented.
#' @param covariance If \code{TRUE}, additionally return (a submatrix of) the generalized inverse of graph Laplacian across the grid
#' @param cores Number of worker processes to use. \code{1} evaluates the grid serially.
#' @param approximation Exploratory approximation used during grid evaluation.
#'   \code{"none"} evaluates the full focal set at each grid point.
#'   \code{"landmark"} evaluates a landmark subset selected once and reused
#'   across the whole grid.
#'   \code{"coarse_raster"} keeps the full focal set but rebuilds the graph on
#'   an aggregated raster for safer coarse screening.
#' @param approximation_control Optional named list controlling the landmark
#'   or coarse-raster approximation. For \code{"landmark"}, see
#'   \code{\link{terradish}} for the supported entries. For
#'   \code{"coarse_raster"}, supported entries include \code{factor},
#'   \code{aggregate_fun}, and \code{directions}. Coarse-raster screening
#'   requires that \code{data} retain its original raster stack.
#'
#' @return An object of class \code{terradish_grid} with components
#'   \code{theta}, \code{loglik}, \code{phi}, and, when requested,
#'   \code{covariance}.
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
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' fit_mlpe <- terradish(melip.Fst ~ altitude + forestcover, data = surface,
#'                    terradish::loglinear_conductance, terradish::mlpe,
#'                    control = NewtonRaphsonControl(maxit = 5, verbose = FALSE))
#'
#' theta <- as.matrix(expand.grid(forestcover = seq(-0.5, 0.5, length.out = 3),
#'                                altitude = seq(-0.5, 0.5, length.out = 3)))
#'
#' grid <- terradish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
#'                     terradish::loglinear_conductance, terradish::mlpe, cores = 1)
#'
#' grid_coarse <- terradish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
#'                            terradish::loglinear_conductance, terradish::mlpe,
#'                            cores = 1,
#'                            approximation = "coarse_raster",
#'                            approximation_control = list(factor = 2L))
#'
#' cbind(grid$theta, loglik = grid$loglik)
#' coef(fit_mlpe)
#'
#' @export

terradish_grid <- function(theta,
                        formula, 
                        data,
                        conductance_model = loglinear_conductance, 
                        measurement_model = mlpe, 
                        nu = NULL, 
                        nonnegative = TRUE, 
                        conductance = TRUE,
                        covariance  = FALSE,
                        cores = 1L,
                        approximation = c("none", "landmark", "coarse_raster"),
                        approximation_control = NULL)
{
  stopifnot(is.matrix(theta))
  stopifnot(length(cores) == 1, is.numeric(cores), cores >= 1)
  if (!isTRUE(conductance))
    stop("`conductance = FALSE` is not currently supported.", call. = FALSE)
  approximation <- match.arg(approximation)
  conductance_model_factory <- conductance_model
  trm <- terms(formula)
  vars <- as.character(attr(trm, "variables"))[-1]
  response <- attr(trm, "response")
  S <- if (response) get(vars[response], parent.frame())
       else stop("'formula' must have genetic distance matrix on lhs")
  stopifnot(length(vars) > 1)

  term_labels <- attr(trm, "term.labels")
  force(theta)
  force(term_labels)
  force(data)
  force(S)
  force(conductance_model_factory)
  force(measurement_model)
  force(nu)
  .terradish_grid_impl(
    theta = theta,
    term_labels = term_labels,
    data = data,
    S = S,
    conductance_model_factory = conductance_model_factory,
    measurement_model = measurement_model,
    nu = nu,
    nonnegative = nonnegative,
    covariance = covariance,
    cores = cores,
    approximation = approximation,
    approximation_control = approximation_control
  )
}

#' @rdname terradish_grid
#' @param ... Arguments passed through the deprecated \code{radish_grid()}
#'   compatibility wrapper to \code{\link{terradish_grid}}.
#' @export
radish_grid <- function(...)
{
  .terradish_deprecate("radish_grid", "terradish_grid")
  .terradish_forward_call(match.call(), "terradish_grid")
}

#' Resistance distances from a parameterized conductance surface
#'
#' Calculates resistance distances associated with a parameterized conductance surface across
#' a grid of parameter values.
#'
#' @param theta A matrix of dimension (grid size) x (number of parameters)
#' @param formula A formula with the name of a matrix of observed genetic distances on the lhs, and covariates in the creation of \code{data} on the rhs
#' @param data An object of class \code{terradish_graph} (see
#'   \code{\link{conductance_surface}})
#' @param conductance_model A function of class
#'   \code{terradish_conductance_model_factory} (see
#'   \code{\link{terradish_conductance_model_factory}})
#' @param conductance Retained for backward compatibility. Only
#'   \code{conductance = TRUE} is currently implemented.
#' @param covariance If \code{TRUE}, instead of a matrix of resistance distances, return the associated submatrix of the generalized inverse of graph Laplacian
#' @param cores Number of worker processes to use. \code{1} evaluates the grid serially.
#'
#' @return An object of class \code{terradish_grid}
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
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#'
#' theta <- as.matrix(expand.grid(forestcover = seq(-0.5, 0.5, length.out = 3),
#'                                altitude = seq(-0.5, 0.5, length.out = 3)))
#'
#' distances <- terradish_distance(theta, ~forestcover + altitude, 
#'                              surface, terradish::loglinear_conductance, cores = 1)
#'
#' ibd <- which(theta[,1] == 0 & theta[,2] == 0)
#' round(distances$distance[,,ibd], 4)
#'
#' @export

terradish_distance <- function(theta,
                            formula, 
                            data,
                            conductance_model = loglinear_conductance, 
                            conductance = TRUE,
                            covariance  = FALSE,
                            cores = 1L)
{
  stopifnot(is.matrix(theta))
  stopifnot(length(cores) == 1, is.numeric(cores), cores >= 1)
  if (!isTRUE(conductance))
    stop("`conductance = FALSE` is not currently supported.", call. = FALSE)
  conductance_model_factory <- conductance_model

  # get response, remove lhs from formula
  terms    <- terms(formula)
  is_ibd   <- length(attr(terms, "factors")) == 0

  stopifnot(!is_ibd) #IBD; nothing to do

  formula  <- reformulate(attr(terms, "term.labels"))

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model <- conductance_model_factory(formula, data$x)
  default <- attr(conductance_model, "default")

  stopifnot(ncol(theta) == length(default))

  theta <- .validate_theta_grid(theta, names(default))
  force(theta)
  force(conductance_model)
  force(conductance_model_factory)
  force(data)

  cv <- array(NA, c(length(data$demes), length(data$demes), nrow(theta)))

  worker_state <- list(theta = theta,
                       conductance_model = conductance_model,
                       conductance_model_factory = conductance_model_factory,
                       conductance_model_factory_name = .parallel_fun_name(
                         conductance_model_factory,
                         c("loglinear_conductance", "linear_conductance")
                       ),
                       data = data,
                       covariance = covariance,
                       term_labels = attr(terms(formula), "term.labels"),
                       radish_algorithm = terradish_algorithm,
                       leastsquares = leastsquares,
                       dist_from_cov = dist_from_cov)

  idx <- seq_len(nrow(theta))
  if (as.integer(cores) <= 1L || nrow(theta) <= 1L)
  {
    results <- .terradish_distance_worker(idx, worker_state)
  }
  else
  {
    splits <- split(idx, cut(idx, breaks = min(as.integer(cores), length(idx)), labels = FALSE))
    cl <- makeCluster(length(splits))
    on.exit(stopCluster(cl), add = TRUE)
    worker_libpaths <- .libPaths()
    use_namespace_workers <- .use_namespace_workers()
    worker_state$use_namespace_workers <- use_namespace_workers
    if (use_namespace_workers)
    {
      clusterExport(cl, varlist = c("worker_libpaths"), envir = environment())
      clusterEvalQ(cl, {
        .libPaths(worker_libpaths)
        library(terradish)
        NULL
      })
    }
    else
    {
      worker_files <- c("R/backtracking.R", "R/hager_zhang.R", "R/newton_raphson.R",
                        "R/radish_conductance_model.R", "R/leastsquares.R",
                        "R/radish_subproblem.R", "R/radish_algorithm.R")
      worker_wd <- getwd()
      clusterExport(cl, varlist = c("worker_files", "worker_wd"), envir = environment())
      clusterEvalQ(cl, {
        setwd(worker_wd)
        for (f in worker_files)
          source(f)
        NULL
      })
    }
    worker_state$radish_algorithm <- if (use_namespace_workers) get("terradish_algorithm", envir = asNamespace("terradish")) else terradish_algorithm
    results <- unlist(parLapply(cl, splits, .terradish_distance_worker, state = worker_state), recursive = FALSE)
  }

  for (res in results)
  {
    if (!is.null(res$error))
      next
    cv[, , res$i] <- res$value
  }

  df <- list(theta = data.frame(theta))
  if (covariance) df$covariance <- cv else df$distance <- cv

  class(df) <- c("terradish_grid", "radish_grid")
  df
}

#' @rdname terradish_distance
#' @param ... Arguments passed through the deprecated
#'   \code{radish_distance()} compatibility wrapper to
#'   \code{\link{terradish_distance}}.
#' @export
radish_distance <- function(...)
{
  .terradish_deprecate("radish_distance", "terradish_distance")
  .terradish_forward_call(match.call(), "terradish_distance")
}

# Internal compatibility aliases retained for housekeeping-sized refactors.
.radish_grid_worker <- .terradish_grid_worker
.radish_distance_worker <- .terradish_distance_worker
.radish_grid_approximation <- .terradish_grid_approximation
.radish_grid_parallel_chunk <- .terradish_grid_parallel_chunk
.radish_grid_parallel_results <- .terradish_grid_parallel_results
.radish_grid_serial_results <- .terradish_grid_serial_results
.radish_grid_impl <- .terradish_grid_impl
