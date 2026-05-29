.terradish_assessment_response <- function(formula, env)
{
  terms_obj <- terms(formula)
  response <- attr(terms_obj, "response")
  if (!response)
    stop("`formula` must have a genetic distance matrix on the left-hand side",
         call. = FALSE)

  vars <- as.character(attr(terms_obj, "variables"))[-1L]
  response_name <- vars[[response]]
  if (!exists(response_name, envir = env, inherits = TRUE))
    stop("Could not find response object `", response_name, "` for assessment",
         call. = FALSE)

  list(name = response_name,
       value = get(response_name, envir = env, inherits = TRUE))
}

.terradish_assessment_rhs_formula <- function(formula)
{
  terms_obj <- terms(formula)
  vars <- as.character(attr(terms_obj, "variables"))[-1L]
  if (length(vars) == 1L)
    return(formula(~1))

  reformulate(attr(terms_obj, "term.labels"))
}

.terradish_assessment_profile <- function(formula, data, conductance_model,
                                          theta = NULL)
{
  rhs_formula <- .terradish_assessment_rhs_formula(formula)
  model <- conductance_model(rhs_formula, data$x)
  default <- attr(model, "default", exact = TRUE)
  if (is.null(default))
    stop("Conductance model does not define default parameters", call. = FALSE)

  theta_internal <- if (is.null(theta))
    default
  else
  {
    if (length(theta) != length(default))
      stop("`theta` must have one value per conductance parameter",
           call. = FALSE)
    names(theta) <- names(.conductance_model_parameter_scale(model))
    .conductance_model_to_internal(theta, model)
  }

  can_coarse <- !is.null(data$stack)
  if (can_coarse)
  {
    can_coarse <- tryCatch({
      .validate_multiscale_covariates(
        data$stack,
        caller = "`terradish_assess_settings()`"
      )
      TRUE
    }, error = function(e) FALSE)
  }

  smooth_info <- attr(model, "smooth_loglinear_info", exact = TRUE)
  smooth_specs <- smooth_info$smooth_specs
  smooth_basis_columns <- if (length(smooth_specs))
    sum(vapply(smooth_specs, function(spec) length(spec$columns), integer(1)))
  else
    0L

  list(
    n_vertices = nrow(data$x),
    n_edges = ncol(data$adj),
    n_focal = length(data$demes),
    n_rhs = ncol(.graph_rhs(data, nrow(data$x))),
    n_parameters = length(default),
    parameter_names = names(default),
    has_stack = !is.null(data$stack),
    can_coarse_raster = isTRUE(can_coarse),
    smooth_conductance = isTRUE(attr(model, "smooth_loglinear", exact = TRUE)),
    n_smooth_terms = length(smooth_specs),
    n_smooth_basis_columns = as.integer(smooth_basis_columns),
    rhs_formula = rhs_formula,
    conductance_model = model,
    theta_internal = theta_internal
  )
}

.terradish_assessment_fit_probe <- function(label,
                                            formula,
                                            response_name,
                                            response,
                                            data,
                                            conductance_model,
                                            measurement_model,
                                            nu,
                                            theta,
                                            nonnegative,
                                            optimizer,
                                            ls_control,
                                            probe_maxit,
                                            solver,
                                            solver_control,
                                            approximation = "none",
                                            approximation_control = NULL,
                                            cores = 1L)
{
  assign(response_name, response, envir = environment())
  control <- NewtonRaphsonControl(
    maxit = probe_maxit,
    verbose = FALSE,
    ctol = 1e-6,
    ftol = 1e-6,
    ls.control = ls_control
  )

  gc()
  elapsed_start <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    suppressWarnings(
      terradish(
        formula,
        data = data,
        conductance_model = conductance_model,
        measurement_model = measurement_model,
        nu = nu,
        theta = theta,
        leverage = FALSE,
        nonnegative = nonnegative,
        optimizer = optimizer,
        control = control,
        cores = cores,
        solver = solver,
        solver_control = solver_control,
        approximation = approximation,
        approximation_control = approximation_control
      )
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[["elapsed"]] - elapsed_start

  if (inherits(fit, "error"))
  {
    return(data.frame(
      label = label,
      optimizer = optimizer,
      line_search = if (inherits(ls_control, "terradish_armijo_control")) "armijo" else "hager_zhang",
      approximation = approximation,
      elapsed = elapsed,
      loglik = NA_real_,
      steps = NA_integer_,
      calls = NA_integer_,
      objective_evaluations = NA_real_,
      gradient_evaluations = NA_real_,
      hessian_evaluations = NA_real_,
      line_search_trials = NA_real_,
      line_search_objective_only_trials = NA_real_,
      solver_setups = NA_real_,
      solver_solves = NA_real_,
      solver_setup_time = NA_real_,
      solver_solve_time = NA_real_,
      status = paste("ERROR:", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  diagnostics <- fit$diagnostics
  data.frame(
    label = label,
    optimizer = optimizer,
    line_search = if (inherits(ls_control, "terradish_armijo_control")) "armijo" else "hager_zhang",
    approximation = approximation,
    elapsed = elapsed,
    loglik = fit$loglik,
    steps = unname(as.integer(fit$cost[["newton_steps"]])),
    calls = unname(as.integer(fit$cost[["function_calls"]])),
    objective_evaluations = diagnostics$objective_evaluations,
    gradient_evaluations = diagnostics$gradient_evaluations,
    hessian_evaluations = diagnostics$hessian_evaluations,
    line_search_trials = diagnostics$line_search_trials,
    line_search_objective_only_trials = diagnostics$line_search_objective_only_trials,
    solver_setups = diagnostics$solver_setups,
    solver_solves = diagnostics$solver_solves,
    solver_setup_time = diagnostics$solver_setup_time,
    solver_solve_time = diagnostics$solver_solve_time,
    status = "OK",
    stringsAsFactors = FALSE
  )
}

.terradish_assessment_candidate_spec <- function(candidate)
{
  switch(candidate,
         newton_hager_zhang = list(optimizer = "newton",
                                   ls_control = HagerZhangControl()),
         bfgs_hager_zhang = list(optimizer = "bfgs",
                                  ls_control = HagerZhangControl()),
         bfgs_armijo = list(optimizer = "bfgs",
                             ls_control = ArmijoControl()),
         stop("Unknown optimizer candidate: ", candidate, call. = FALSE))
}

.terradish_assessment_pick_fit <- function(fit_benchmark,
                                           fallback_optimizer,
                                           loglik_tolerance)
{
  if (is.null(fit_benchmark) || !nrow(fit_benchmark))
    return(list(optimizer = fallback_optimizer,
                line_search = "hager_zhang",
                reason = "profile heuristic"))

  ok <- fit_benchmark[fit_benchmark$status == "OK", , drop = FALSE]
  if (!nrow(ok))
    return(list(optimizer = fallback_optimizer,
                line_search = "hager_zhang",
                reason = "all optimizer probes failed"))

  best_loglik <- max(ok$loglik, na.rm = TRUE)
  acceptable <- ok[ok$loglik >= best_loglik - loglik_tolerance, , drop = FALSE]
  if (!nrow(acceptable))
    acceptable <- ok
  best <- acceptable[which.min(acceptable$elapsed), , drop = FALSE]

  list(optimizer = best$optimizer[[1]],
       line_search = best$line_search[[1]],
       reason = paste0("fastest acceptable probe: ", best$label[[1]]))
}

.terradish_assessment_pick_solver <- function(profile, solver_benchmark)
{
  if (is.null(solver_benchmark) || !nrow(solver_benchmark))
  {
    solver <- if (profile$n_vertices >= 1500000L && profile$n_rhs <= 64L)
      "auto"
    else
      "direct"
    return(list(solver = solver,
                solver_control = NULL,
                reason = "profile heuristic"))
  }

  best <- solver_benchmark[which.min(solver_benchmark$total_time), , drop = FALSE]
  list(
    solver = "direct",
    solver_control = list(
      factorization = best$factorization[[1]],
      solve_backend = best$solve_backend[[1]]
    ),
    reason = paste0("fastest direct solver probe: ",
                    best$factorization[[1]], " + ",
                    best$solve_backend[[1]])
  )
}

.terradish_assessment_control <- function(line_search)
{
  ls_control <- if (identical(line_search, "armijo"))
    ArmijoControl()
  else
    HagerZhangControl()

  NewtonRaphsonControl(
    verbose = FALSE,
    ctol = 1e-6,
    ftol = 1e-6,
    ls.control = ls_control
  )
}

.terradish_assessment_line_search_name <- function(control)
{
  if (inherits(control, "terradish_armijo_control"))
    "armijo"
  else
    "hager_zhang"
}

.terradish_assessment_default_settings <- function(profile,
                                                   conductance_model_factory)
{
  list(
    optimizer = .resolve_terradish_optimizer(
      "auto",
      profile$n_parameters,
      conductance_model_factory = conductance_model_factory
    ),
    control = .terradish_assessment_control("hager_zhang"),
    solver = "direct",
    solver_control = NULL,
    approximation = "none",
    approximation_control = NULL
  )
}

.terradish_assessment_compare_settings <- function(defaults, recommended)
{
  default_line_search <- .terradish_assessment_line_search_name(
    defaults$control$ls.control
  )
  recommended_line_search <- .terradish_assessment_line_search_name(
    recommended$control$ls.control
  )

  changes <- character()
  if (!identical(defaults$optimizer, recommended$optimizer))
    changes <- c(
      changes,
      paste0("optimizer `", defaults$optimizer, "` -> `",
             recommended$optimizer, "`")
    )
  if (!identical(default_line_search, recommended_line_search))
    changes <- c(
      changes,
      paste0("line search `", default_line_search, "` -> `",
             recommended_line_search, "`")
    )
  if (!identical(defaults$solver, recommended$solver))
    changes <- c(
      changes,
      paste0("solver `", defaults$solver, "` -> `",
             recommended$solver, "`")
    )
  if (!identical(defaults$solver_control, recommended$solver_control))
    changes <- c(changes, "solver_control tuned from terradish defaults")
  if (!identical(defaults$approximation, recommended$approximation))
    changes <- c(
      changes,
      paste0("approximation `", defaults$approximation, "` -> `",
             recommended$approximation, "`")
    )

  differs <- length(changes) > 0L
  summary <- if (differs)
    paste0("Assessment differs from terradish defaults: ",
           paste(changes, collapse = "; "), ".")
  else
    "Assessment matched terradish defaults."

  list(
    default_optimizer = defaults$optimizer,
    recommended_optimizer = recommended$optimizer,
    default_line_search = default_line_search,
    recommended_line_search = recommended_line_search,
    default_solver = defaults$solver,
    recommended_solver = recommended$solver,
    default_approximation = defaults$approximation,
    recommended_approximation = recommended$approximation,
    differs_from_defaults = differs,
    changes = changes,
    summary = summary
  )
}

#' Assess terradish speed settings for a data set
#'
#' Profiles a \code{terradish} graph and, optionally, runs short probe
#' benchmarks to recommend optimizer, line-search, solver, and approximation
#' settings for a fitted model.  Use this function before a large or slow fit
#' to identify the fastest computational configuration for your data.
#'
#' @param formula Model formula passed to \code{\link{terradish}}.
#' @param data A \code{terradish_graph} object returned by
#'   \code{\link{conductance_surface}}.
#' @param conductance_model Conductance-model factory (e.g.
#'   \code{\link{loglinear_conductance}}).
#' @param measurement_model Measurement model used during optimizer probes.
#'   Should match the model you plan to use for the final fit.
#' @param nu Effective Wishart degrees of freedom passed to measurement models
#'   that require it (\code{\link{generalized_wishart}},
#'   \code{\link{wishart_covariance}}).  For biallelic SNPs this is usually the
#'   number of retained SNPs; for microsatellites use approximately
#'   \eqn{\sum_l (K_l - 1)} where \eqn{K_l} is the number of observed alleles at
#'   locus \eqn{l}.  Ignored by \code{\link{mlpe}} and
#'   \code{\link{leastsquares}}.
#' @param theta Optional starting conductance parameter values.  If
#'   \code{NULL}, the default all-zero start is used.
#' @param nonnegative Should regression measurement models constrain the IBR
#'   slope to be nonnegative during probes?
#' @param optimizer_probe Logical.  Run short optimization probes to compare
#'   Newton-Raphson vs. BFGS and line-search algorithms?
#' @param solver_probe Logical.  Benchmark direct Cholesky solver settings
#'   (factorization style and solve backend) via
#'   \code{\link{terradish_solver_benchmark}}?
#' @param coarse_probe Logical.  Include a coarse-raster warm-start probe?
#'   Only applicable when \code{data} retains its raster stack
#'   (\code{saveStack = TRUE}).
#' @param optimizer_candidates Character vector of optimizer configurations to
#'   probe.  Supported entries:
#'   \code{"newton_hager_zhang"}, \code{"bfgs_hager_zhang"},
#'   \code{"bfgs_armijo"}.
#' @param probe_maxit Maximum optimizer steps per probe.  Keep small (default
#'   \code{2}); the assessment compares early iteration cost, not convergence.
#' @param factorization Character vector of Cholesky factorization styles to
#'   benchmark when \code{solver_probe = TRUE}.  Options:
#'   \code{"simplicial_ldl"}, \code{"simplicial_ll"}, \code{"supernodal_ll"}.
#' @param solve_backend Character vector of solve backends to benchmark when
#'   \code{solver_probe = TRUE}.  Options:
#'   \code{"matrix"}, \code{"cholmod_cpp"}, \code{"cholmod_cpp_cached"}.
#' @param coarse_factor Raster aggregation factor for the coarse-raster probe
#'   (passed to the \code{factor} entry of \code{approximation_control}).
#' @param coarse_refine_maxit Maximum full-resolution refinement steps for the
#'   coarse-raster probe.
#' @param loglik_tolerance Maximum log-likelihood difference below the best
#'   probe that is still considered acceptable when choosing the fastest option.
#' @param max_seconds Soft wall-time budget in seconds.  Each probe that starts
#'   before the budget is exceeded is allowed to finish; later probes are
#'   skipped.
#' @param cores Number of parallel worker processes passed to probes.
#' @param verbose Logical.  Print progress as probes run?
#'
#' @details
#' \code{terradish_assess_settings} is an advisory tool, not a model-selection
#' procedure.  It does not change package defaults and its output does not
#' affect parameter estimation.  Use it to choose \emph{computational}
#' settings; then fit and interpret your model as usual.
#'
#' Short probes work well for comparing early-iteration cost of different
#' optimizers.  They are most informative for large graphs (> 100,000 vertices)
#' where a single function evaluation is slow.  For small graphs the default
#' Newton-Raphson + Hager-Zhang settings typically converge in seconds, and
#' assessment adds more overhead than it saves.
#'
#' The returned \code{$recommended} list is intentionally a plain named list of
#' arguments (\code{optimizer}, \code{control}, \code{solver},
#' \code{solver_control}, \code{approximation},
#' \code{approximation_control}).  You can inspect, edit, or pass it directly:
#' \preformatted{
#' fit <- terradish(formula, data,
#'                  optimizer = assessment$recommended$optimizer,
#'                  control   = assessment$recommended$control,
#'                  solver    = assessment$recommended$solver)
#' }
#'
#' @return An object of class \code{"terradish_setting_assessment"} with
#'   components:
#' \describe{
#'   \item{\code{profile}}{A list of graph statistics (number of vertices,
#'     edges, focal sites, parameters, and, for
#'     \code{\link{smooth_loglinear_conductance}}, spline-expansion counts.}
#'   \item{\code{optimizer_benchmark}}{A data frame of optimizer probe results
#'     (one row per candidate), or \code{NULL} if \code{optimizer_probe = FALSE}.}
#'   \item{\code{solver_benchmark}}{A data frame of solver probe results, or
#'     \code{NULL} if \code{solver_probe = FALSE}.}
#'   \item{\code{defaults}}{The terradish defaults for the profiled model,
#'     after resolving \code{optimizer = "auto"} on the expanded conductance
#'     parameter count.}
#'   \item{\code{comparison}}{A compact comparison between the recommended
#'     settings and the terradish defaults, including a one-line summary.}
#'   \item{\code{recommended}}{A named list of recommended settings:
#'     \code{optimizer}, \code{control}, \code{solver}, \code{solver_control},
#'     \code{approximation}, \code{approximation_control}.}
#'   \item{\code{notes}}{Character vector of advisory messages.}
#' }
#' A \code{print} method summarizes the profile and recommendations.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' covariates <- c(melip.altitude, melip.forestcover)
#' names(covariates) <- c("altitude", "forestcover")
#' covariates <- scale_covariates(covariates)
#' surface <- conductance_surface(covariates, melip.coords,
#'                                directions = 8, saveStack = TRUE)
#'
#' assessment <- terradish_assess_settings(
#'   melip.Fst ~ forestcover + altitude,
#'   data = surface,
#'   measurement_model = mlpe,
#'   probe_maxit = 2
#' )
#' assessment
#'
#' fit <- terradish(
#'   melip.Fst ~ forestcover + altitude,
#'   data = surface,
#'   conductance_model = loglinear_conductance,
#'   measurement_model = mlpe,
#'   optimizer = assessment$recommended$optimizer,
#'   control = assessment$recommended$control,
#'   solver = assessment$recommended$solver,
#'   solver_control = assessment$recommended$solver_control
#' )
#'
#' smooth_assessment <- terradish_assess_settings(
#'   melip.Fst ~ forestcover + s(altitude, df = 3),
#'   data = surface,
#'   conductance_model = smooth_loglinear_conductance,
#'   measurement_model = mlpe,
#'   probe_maxit = 1
#' )
#' smooth_assessment$comparison$summary
#' }
#'
#' @export
terradish_assess_settings <- function(formula,
                                      data,
                                      conductance_model = loglinear_conductance,
                                      measurement_model = mlpe,
                                      nu = NULL,
                                      theta = NULL,
                                      nonnegative = TRUE,
                                      optimizer_probe = TRUE,
                                      solver_probe = TRUE,
                                      coarse_probe = FALSE,
                                      optimizer_candidates = c("newton_hager_zhang",
                                                               "bfgs_hager_zhang",
                                                               "bfgs_armijo"),
                                      probe_maxit = 2L,
                                      factorization = c("simplicial_ldl",
                                                        "simplicial_ll",
                                                        "supernodal_ll"),
                                      solve_backend = c("matrix",
                                                        "cholmod_cpp",
                                                        "cholmod_cpp_cached"),
                                      coarse_factor = 2L,
                                      coarse_refine_maxit = 1L,
                                      loglik_tolerance = 1e-3,
                                      max_seconds = 120,
                                      cores = 1L,
                                      verbose = TRUE)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(inherits(conductance_model, c("terradish_conductance_model_factory",
                                          "radish_conductance_model_factory")))
  stopifnot(inherits(measurement_model, c("terradish_measurement_model",
                                          "radish_measurement_model")))

  probe_maxit <- as.integer(probe_maxit)[1]
  if (is.na(probe_maxit) || probe_maxit < 1L)
    stop("`probe_maxit` must be a positive integer", call. = FALSE)
  if (!is.finite(loglik_tolerance) || loglik_tolerance < 0)
    stop("`loglik_tolerance` must be a nonnegative finite value", call. = FALSE)

  caller <- parent.frame()
  response <- .terradish_assessment_response(formula, caller)
  profile <- .terradish_assessment_profile(formula, data, conductance_model,
                                           theta = theta)
  notes <- character()
  started <- proc.time()[["elapsed"]]
  over_budget <- function()
    is.finite(max_seconds) &&
      (proc.time()[["elapsed"]] - started) > max_seconds

  fallback_optimizer <- .resolve_terradish_optimizer(
    "auto",
    profile$n_parameters,
    conductance_model_factory = conductance_model
  )
  defaults <- .terradish_assessment_default_settings(
    profile,
    conductance_model_factory = conductance_model
  )

  if (isTRUE(verbose))
    message("Profiling graph: ", profile$n_vertices, " vertices, ",
            profile$n_focal, " focal sites, ",
            profile$n_parameters, " conductance parameters")

  solver_benchmark <- NULL
  if (isTRUE(solver_probe) && !over_budget())
  {
    if (isTRUE(verbose))
      message("Benchmarking direct solver settings")
    conductance <- profile$conductance_model(profile$theta_internal)$conductance
    solver_benchmark <- terradish_solver_benchmark(
      data,
      conductance = conductance,
      factorization = factorization,
      solve_backend = solve_backend,
      n_replicates = 1L
    )
    solver_benchmark <- solver_benchmark[order(solver_benchmark$total_time), ,
                                         drop = FALSE]
  }
  else if (isTRUE(solver_probe))
    notes <- c(notes, "Skipped solver probes because the soft time budget was exceeded.")

  fit_benchmark <- NULL
  if (isTRUE(optimizer_probe) && !over_budget())
  {
    optimizer_candidates <- match.arg(
      optimizer_candidates,
      c("newton_hager_zhang", "bfgs_hager_zhang", "bfgs_armijo"),
      several.ok = TRUE
    )
    if (isTRUE(verbose))
      message("Running short optimizer probes")
    fit_rows <- vector("list", length(optimizer_candidates))
    for (i in seq_along(optimizer_candidates))
    {
      candidate <- optimizer_candidates[[i]]
      spec <- .terradish_assessment_candidate_spec(candidate)
      fit_rows[[i]] <- .terradish_assessment_fit_probe(
        label = candidate,
        formula = formula,
        response_name = response$name,
        response = response$value,
        data = data,
        conductance_model = conductance_model,
        measurement_model = measurement_model,
        nu = nu,
        theta = theta,
        nonnegative = nonnegative,
        optimizer = spec$optimizer,
        ls_control = spec$ls_control,
        probe_maxit = probe_maxit,
        solver = "direct",
        solver_control = .terradish_assessment_pick_solver(profile, solver_benchmark)$solver_control,
        cores = cores
      )
      if (over_budget())
      {
        notes <- c(notes, "Stopped optimizer probes after the soft time budget was exceeded.")
        fit_rows <- fit_rows[seq_len(i)]
        break
      }
    }
    fit_benchmark <- do.call(rbind, fit_rows)
  }
  else if (isTRUE(optimizer_probe))
    notes <- c(notes, "Skipped optimizer probes because the soft time budget was exceeded.")

  coarse_benchmark <- NULL
  if (isTRUE(coarse_probe) && !over_budget())
  {
    if (!isTRUE(profile$can_coarse_raster))
    {
      notes <- c(notes, "Skipped coarse-raster probe because data does not retain a continuous raster stack.")
    }
    else
    {
      if (isTRUE(verbose))
        message("Running short coarse-raster warm-start probe")
      fit_pick <- .terradish_assessment_pick_fit(
        fit_benchmark,
        fallback_optimizer,
        loglik_tolerance = loglik_tolerance
      )
      coarse_benchmark <- .terradish_assessment_fit_probe(
        label = paste0("coarse_factor_", coarse_factor),
        formula = formula,
        response_name = response$name,
        response = response$value,
        data = data,
        conductance_model = conductance_model,
        measurement_model = measurement_model,
        nu = nu,
        theta = theta,
        nonnegative = nonnegative,
        optimizer = fit_pick$optimizer,
        ls_control = if (identical(fit_pick$line_search, "armijo"))
          ArmijoControl()
        else
          HagerZhangControl(),
        probe_maxit = probe_maxit,
        solver = "direct",
        solver_control = .terradish_assessment_pick_solver(profile, solver_benchmark)$solver_control,
        approximation = "coarse_raster",
        approximation_control = list(
          factor = as.integer(coarse_factor)[1],
          exact_refine = TRUE,
          refine_control = NewtonRaphsonControl(
            maxit = as.integer(coarse_refine_maxit)[1],
            verbose = FALSE
          )
        ),
        cores = cores
      )
    }
  }

  fit_pick <- .terradish_assessment_pick_fit(
    fit_benchmark,
    fallback_optimizer,
    loglik_tolerance = loglik_tolerance
  )
  solver_pick <- .terradish_assessment_pick_solver(profile, solver_benchmark)

  approximation <- "none"
  approximation_control <- NULL
  if (!is.null(coarse_benchmark) &&
      nrow(coarse_benchmark) &&
      identical(coarse_benchmark$status[[1]], "OK") &&
      !is.null(fit_benchmark) &&
      any(fit_benchmark$status == "OK"))
  {
    exact_ok <- fit_benchmark[fit_benchmark$status == "OK", , drop = FALSE]
    best_exact_loglik <- max(exact_ok$loglik, na.rm = TRUE)
    best_exact_elapsed <- min(exact_ok$elapsed, na.rm = TRUE)
    if (coarse_benchmark$elapsed[[1]] < best_exact_elapsed &&
        coarse_benchmark$loglik[[1]] >= best_exact_loglik - loglik_tolerance)
    {
      approximation <- "coarse_raster"
      approximation_control <- list(
        factor = as.integer(coarse_factor)[1],
        exact_refine = TRUE,
        refine_control = NewtonRaphsonControl(
          maxit = as.integer(coarse_refine_maxit)[1],
          verbose = FALSE
        )
      )
      notes <- c(notes, "Coarse-raster probe was faster and within the log-likelihood tolerance.")
    }
    else
      notes <- c(notes, "Coarse-raster probe was not recommended over the exact probe.")
  }

  recommended <- list(
    optimizer = fit_pick$optimizer,
    control = .terradish_assessment_control(fit_pick$line_search),
    solver = solver_pick$solver,
    solver_control = solver_pick$solver_control,
    approximation = approximation,
    approximation_control = approximation_control
  )
  reasons <- list(
    optimizer = fit_pick$reason,
    solver = solver_pick$reason,
    approximation = if (identical(approximation, "none"))
      "exact full-resolution fit"
    else
      "coarse probe accepted"
  )
  comparison <- .terradish_assessment_compare_settings(defaults, recommended)

  if (isTRUE(profile$smooth_conductance))
  {
    notes <- c(
      notes,
      paste0(
        "Smooth conductance formula expanded to ",
        profile$n_parameters,
        " conductance parameters across ",
        profile$n_smooth_terms,
        " smooth term",
        if (identical(profile$n_smooth_terms, 1L)) "" else "s",
        " (",
        profile$n_smooth_basis_columns,
        " spline basis columns)."
      )
    )

    if (!is.null(fit_benchmark) && nrow(fit_benchmark))
    {
      failed <- fit_benchmark[fit_benchmark$status != "OK", , drop = FALSE]
      if (nrow(failed))
      {
        notes <- c(
          notes,
          paste0(
            "Some smooth-model probes hit optimizer or conductance guards: ",
            paste(paste0(failed$label, " [", failed$status, "]"),
                  collapse = "; "),
            "."
          )
        )
      }
      else
        notes <- c(
          notes,
          "Short smooth-model probes completed without additional optimizer or conductance guard failures."
        )
    }
  }

  out <- list(
    call = match.call(),
    profile = profile[setdiff(names(profile),
                              c("conductance_model", "theta_internal"))],
    defaults = defaults,
    comparison = comparison,
    recommended = recommended,
    reasons = reasons,
    benchmarks = list(
      solver = solver_benchmark,
      optimizer = fit_benchmark,
      coarse = coarse_benchmark
    ),
    notes = notes
  )
  class(out) <- "terradish_setting_assessment"
  out
}

#' @method print terradish_setting_assessment
#' @export
print.terradish_setting_assessment <- function(x,
                                               digits = max(3L, getOption("digits") - 3L),
                                               ...)
{
  cat("terradish setting assessment\n")
  cat("Graph:",
      x$profile$n_vertices, "vertices,",
      x$profile$n_edges, "edges,",
      x$profile$n_focal, "focal sites\n")
  cat("Conductance parameters:", x$profile$n_parameters, "\n\n")

  cat("Recommended settings:\n")
  cat("  optimizer: ", x$recommended$optimizer,
      " (", x$reasons$optimizer, ")\n", sep = "")
  line_search <- if (inherits(x$recommended$control$ls.control,
                              "terradish_armijo_control"))
    "ArmijoControl"
  else
    "HagerZhangControl"
  cat("  line search: ", line_search, "\n", sep = "")
  cat("  solver: ", x$recommended$solver,
      " (", x$reasons$solver, ")\n", sep = "")
  if (!is.null(x$recommended$solver_control))
  {
    cat("  solver_control:\n")
    print(x$recommended$solver_control)
  }
  cat("  approximation: ", x$recommended$approximation,
      " (", x$reasons$approximation, ")\n\n", sep = "")
  cat("Defaults comparison: ", x$comparison$summary, "\n\n", sep = "")

  if (!is.null(x$benchmarks$optimizer))
  {
    cat("Optimizer probes:\n")
    print(format(x$benchmarks$optimizer[, c("label", "elapsed", "loglik",
                                            "gradient_evaluations",
                                            "hessian_evaluations",
                                            "line_search_trials",
                                            "status"),
                                      drop = FALSE],
                 digits = digits),
          row.names = FALSE)
    cat("\n")
  }

  if (!is.null(x$benchmarks$solver))
  {
    cat("Fastest direct solver probes:\n")
    print(format(utils::head(x$benchmarks$solver[, c("factorization",
                                                     "solve_backend",
                                                     "setup_time",
                                                     "solve_time",
                                                     "total_time"),
                                               drop = FALSE], 3L),
                 digits = digits),
          row.names = FALSE)
    cat("\n")
  }

  if (length(x$notes))
  {
    cat("Notes:\n")
    for (note in x$notes)
      cat("  - ", note, "\n", sep = "")
  }

  invisible(x)
}
