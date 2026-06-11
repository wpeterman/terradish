setRefClass("FunctionCall", fields = list(count = "integer"))

.terradish_new_diagnostics <- function()
{
  diagnostics <- new.env(parent = emptyenv())
  for (name in c("algorithm_calls",
                 "objective_evaluations",
                 "gradient_evaluations",
                 "hessian_evaluations",
                 "partial_evaluations",
                 "line_search_trials",
                 "line_search_objective_only_trials",
                 "line_search_gradient_trials",
                 "solver_setups",
                 "solver_solves"))
    diagnostics[[name]] <- 0L
  diagnostics$solver_setup_time <- 0
  diagnostics$solver_solve_time <- 0
  diagnostics
}

.terradish_increment_diagnostic <- function(diagnostics, name, value = 1L)
{
  if (is.null(diagnostics) || !is.environment(diagnostics))
    return(invisible(NULL))
  current <- diagnostics[[name]]
  if (is.null(current))
    current <- if (is.integer(value)) 0L else 0
  diagnostics[[name]] <- current + value
  invisible(diagnostics[[name]])
}

.terradish_record_line_search_trial <- function(control,
                                                objective_only = FALSE,
                                                gradient = FALSE)
{
  diagnostics <- control$diagnostics
  .terradish_increment_diagnostic(diagnostics, "line_search_trials")
  if (isTRUE(objective_only))
    .terradish_increment_diagnostic(diagnostics, "line_search_objective_only_trials")
  if (isTRUE(gradient))
    .terradish_increment_diagnostic(diagnostics, "line_search_gradient_trials")
  invisible(NULL)
}

.terradish_record_algorithm_diagnostics <- function(diagnostics,
                                                    fit,
                                                    objective = TRUE,
                                                    gradient = TRUE,
                                                    hessian = TRUE,
                                                    partial = FALSE)
{
  .terradish_increment_diagnostic(diagnostics, "algorithm_calls")
  if (isTRUE(objective))
    .terradish_increment_diagnostic(diagnostics, "objective_evaluations")
  if (isTRUE(gradient))
    .terradish_increment_diagnostic(diagnostics, "gradient_evaluations")
  if (isTRUE(hessian))
    .terradish_increment_diagnostic(diagnostics, "hessian_evaluations")
  if (isTRUE(partial))
    .terradish_increment_diagnostic(diagnostics, "partial_evaluations")

  solver_diagnostics <- fit$algorithm_diagnostics
  if (is.null(solver_diagnostics))
    return(invisible(NULL))

  .terradish_increment_diagnostic(diagnostics, "solver_setups",
                                  solver_diagnostics$solver_setups)
  .terradish_increment_diagnostic(diagnostics, "solver_solves",
                                  solver_diagnostics$solver_solves)
  .terradish_increment_diagnostic(diagnostics, "solver_setup_time",
                                  solver_diagnostics$solver_setup_time)
  .terradish_increment_diagnostic(diagnostics, "solver_solve_time",
                                  solver_diagnostics$solver_solve_time)
  invisible(NULL)
}

.terradish_diagnostics_snapshot <- function(diagnostics)
{
  if (is.null(diagnostics) || !is.environment(diagnostics))
    return(NULL)

  names <- c("algorithm_calls",
             "objective_evaluations",
             "gradient_evaluations",
             "hessian_evaluations",
             "partial_evaluations",
             "line_search_trials",
             "line_search_objective_only_trials",
             "line_search_gradient_trials",
             "solver_setups",
             "solver_solves",
             "solver_setup_time",
             "solver_solve_time")
  out <- lapply(names, function(name) diagnostics[[name]])
  names(out) <- names
  out
}

.safe_hessian_inverse <- function(x)
{
  tryCatch(solve(x), error = function(e) ginv(x))
}

.conditional_phi_table <- function(phi, phi_hessian, conf.level = 0.95)
{
  if (is.null(phi) || is.null(phi_hessian))
    return(NULL)

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("`conf.level` must be a single number between 0 and 1")

  phi <- c(phi)
  vcov <- .safe_hessian_inverse(phi_hessian)
  se <- sqrt(pmax(diag(vcov), 0))
  z <- qnorm((1 + conf.level) / 2)
  pct <- format(100 * conf.level, trim = TRUE, scientific = FALSE)

  table <- cbind("Estimate" = phi,
                 "Std. Error" = se,
                 phi - z * se,
                 phi + z * se)
  colnames(table) <- c("Estimate", "Std. Error",
                       paste0("Lower ", pct, "%"),
                       paste0("Upper ", pct, "%"))
  rownames(table) <- names(phi)

  list(table = table, vcov = vcov)
}

.terradish_amg_control_defaults <- function()
{
  list(
    adaptive = TRUE,
    tol_early = 1e-4,
    tol_mid = 1e-6,
    tol_final = 1e-8,
    maxit_early = 100L,
    maxit_mid = 250L,
    maxit_final = 400L,
    warmup_evals = 3L
  )
}

.terradish_solver_control_for_phase <- function(solver, solver_control = NULL, eval_count = 0L, final = FALSE)
{
  if (!solver %in% c("amg", "auto"))
    return(solver_control)

  control <- if (is.null(solver_control)) list() else as.list(solver_control)
  defaults <- .terradish_amg_control_defaults()
  control <- modifyList(defaults, control)

  if (!isTRUE(control$adaptive))
  {
    if (is.null(control$tol) && !is.null(control$tol_final))
      control$tol <- control$tol_final
    if (is.null(control$maxit) && !is.null(control$maxit_final))
      control$maxit <- control$maxit_final
    return(control)
  }

  phase <- if (isTRUE(final))
    "final"
  else if (as.integer(eval_count) <= as.integer(control$warmup_evals))
    "early"
  else
    "mid"

  control$tol <- switch(
    phase,
    early = control$tol_early,
    mid = control$tol_mid,
    final = control$tol_final
  )
  control$maxit <- switch(
    phase,
    early = as.integer(control$maxit_early),
    mid = as.integer(control$maxit_mid),
    final = as.integer(control$maxit_final)
  )
  control$adaptive_phase <- phase
  control
}

.terradish_landmark_control_defaults <- function()
{
  list(
    n_landmarks = NULL,
    fraction = 0.5,
    min_landmarks = 8L,
    max_landmarks = 64L,
    method = "spacefill",
    seed = NULL,
    exact_refine = TRUE
  )
}

.normalize_landmark_control <- function(control, n_focal)
{
  control <- if (is.null(control)) list() else as.list(control)
  defaults <- .terradish_landmark_control_defaults()
  control <- modifyList(defaults, control)

  requested <- control$n_landmarks
  if (is.null(requested))
  {
    requested <- ceiling(n_focal * control$fraction)
    requested <- max(as.integer(control$min_landmarks), as.integer(requested))
    requested <- min(as.integer(control$max_landmarks), requested)
  }

  requested <- min(as.integer(control$max_landmarks), requested, as.integer(n_focal))
  control$n_landmarks <- as.integer(requested)
  control$method <- match.arg(control$method, c("spacefill", "random", "sequential"))
  control
}

.terradish_measurement_model_name <- function(model)
{
  known <- c("leastsquares", "mlpe", "generalized_wishart", "wishart_covariance")
  ns_env <- asNamespace("terradish")
  for (nm in known)
  {
    obj <- get0(nm, envir = ns_env, mode = "function", inherits = FALSE)
    if (!is.null(obj) &&
        identical(formals(model), formals(obj)) &&
        identical(body(model), body(obj)))
      return(nm)
  }

  base_model <- attr(model, "base_model", exact = TRUE)
  if (is.character(base_model) && length(base_model) == 1L)
    return(base_model)

  NULL
}

.deme_coordinates <- function(data)
{
  coords <- data$vertex_coordinates
  if (is.null(coords) || is.null(data$demes))
    return(NULL)
  coords <- as.matrix(coords)
  if (nrow(coords) < max(data$demes))
    return(NULL)
  coords[data$demes, , drop = FALSE]
}

.spacefill_landmarks <- function(coords, k)
{
  n <- nrow(coords)
  if (k >= n)
    return(seq_len(n))

  center <- colMeans(coords, na.rm = TRUE)
  d_center <- rowSums(sweep(coords, 2, center, FUN = "-")^2)
  picked <- integer(k)
  picked[1] <- which.min(d_center)

  min_dist <- rowSums(sweep(coords, 2, coords[picked[1], ], FUN = "-")^2)
  min_dist[picked[1]] <- -Inf

  for (i in 2:k)
  {
    picked[i] <- which.max(min_dist)
    d_new <- rowSums(sweep(coords, 2, coords[picked[i], ], FUN = "-")^2)
    min_dist <- pmin(min_dist, d_new)
    min_dist[picked[seq_len(i)]] <- -Inf
  }

  sort(unique(as.integer(picked)))
}

.landmark_indices <- function(data, n_focal, control)
{
  if (control$n_landmarks >= n_focal)
    return(seq_len(n_focal))

  if (!is.null(control$seed))
  {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (old_seed_exists)
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(control$seed)
  }

  idx <- switch(
    control$method,
    random = sort(sample.int(n_focal, control$n_landmarks)),
    sequential = unique(as.integer(round(seq(1, n_focal, length.out = control$n_landmarks)))),
    spacefill = {
      coords <- .deme_coordinates(data)
      if (is.null(coords))
        sort(sample.int(n_focal, control$n_landmarks))
      else
        .spacefill_landmarks(coords, control$n_landmarks)
    }
  )

  if (length(idx) < control$n_landmarks)
  {
    missing <- setdiff(seq_len(n_focal), idx)
    idx <- sort(c(idx, missing[seq_len(control$n_landmarks - length(idx))]))
  }

  idx
}

.terradish_landmark_subset <- function(data, S, approximation_control = NULL)
{
  n_focal <- nrow(S)
  control <- .normalize_landmark_control(approximation_control, n_focal)
  idx <- .landmark_indices(data, n_focal, control)

  subset_data <- data
  subset_data$demes <- data$demes[idx]
  if (!is.null(data$rhs))
    subset_data$rhs <- data$rhs[, idx, drop = FALSE]

  list(
    data = subset_data,
    S = S[idx, idx, drop = FALSE],
    index = idx,
    control = control,
    used = length(idx) < n_focal
  )
}

.conductance_model_bounds <- function(conductance_model)
{
  default <- attr(conductance_model, "default", exact = TRUE)
  if (is.null(default))
    stop("Conductance model does not define default parameters", call. = FALSE)

  lower <- attr(conductance_model, "lower", exact = TRUE)
  upper <- attr(conductance_model, "upper", exact = TRUE)

  if (is.null(lower))
    lower <- rep(-Inf, length(default))
  if (is.null(upper))
    upper <- rep(Inf, length(default))

  lower <- c(lower)
  upper <- c(upper)
  if (length(lower) != length(default) || length(upper) != length(default))
    stop("Conductance-model bounds must have the same length as the default parameter vector",
         call. = FALSE)

  names(lower) <- names(default)
  names(upper) <- names(default)
  list(lower = lower, upper = upper)
}

.conductance_model_parameter_scale <- function(conductance_model)
{
  default <- attr(conductance_model, "default", exact = TRUE)
  if (is.null(default))
    stop("Conductance model does not define default parameters", call. = FALSE)

  scale <- attr(conductance_model, "parameter_scale", exact = TRUE)
  if (is.null(scale))
    scale <- rep(1, length(default))

  scale <- c(scale)
  if (length(scale) != length(default))
    stop("Conductance-model parameter scales must match the parameter vector length",
         call. = FALSE)
  names(scale) <- names(default)
  scale
}

.conductance_model_to_internal <- function(theta, conductance_model)
{
  theta <- c(theta)
  scale <- .conductance_model_parameter_scale(conductance_model)
  names(theta) <- names(scale)
  theta / scale
}

.conductance_model_to_external <- function(theta, conductance_model)
{
  theta <- c(theta)
  scale <- .conductance_model_parameter_scale(conductance_model)
  names(theta) <- names(scale)
  theta * scale
}

.conductance_model_gradient_to_external <- function(gradient, conductance_model)
{
  gradient <- c(gradient)
  scale <- .conductance_model_parameter_scale(conductance_model)
  names(gradient) <- names(scale)
  gradient / scale
}

.conductance_model_hessian_to_external <- function(hessian, conductance_model)
{
  scale <- .conductance_model_parameter_scale(conductance_model)
  sweep(sweep(hessian, 1, scale, "/"), 2, scale, "/")
}

.conductance_model_vcov_to_internal <- function(vcov, conductance_model)
{
  scale <- .conductance_model_parameter_scale(conductance_model)
  sweep(sweep(vcov, 1, scale, "/"), 2, scale, "/")
}

.externalize_conductance_model <- function(conductance_model)
{
  parameter_scale <- .conductance_model_parameter_scale(conductance_model)

  external_model <- function(theta)
  {
    theta_internal <- c(theta) / parameter_scale
    names(theta_internal) <- names(parameter_scale)
    out <- conductance_model(theta_internal)

    out$df__dtheta <- function(k)
      conductance_model(theta_internal)$df__dtheta(k) / parameter_scale[k]
    out$df__dtheta_matrix <- sweep(out$df__dtheta_matrix, 2, parameter_scale, "/")
    out$d2f__dtheta_dtheta <- function(k, l)
      conductance_model(theta_internal)$d2f__dtheta_dtheta(k, l) /
      (parameter_scale[k] * parameter_scale[l])

    confint_internal <- out$confint
    out$confint <- function(theta, vcov, quantile = 0.95,
                            scale = c("conductance", "linpred"))
    {
      scale_arg <- match.arg(scale)
      theta_internal_local <- c(theta) / parameter_scale
      names(theta_internal_local) <- names(parameter_scale)
      vcov_internal <- .conductance_model_vcov_to_internal(vcov, conductance_model)
      confint_internal(theta_internal_local, vcov_internal,
                       quantile = quantile, scale = scale_arg)
    }

    out
  }

  class(external_model) <- class(conductance_model)
  attrs <- attributes(conductance_model)
  if (!is.null(attr(conductance_model, "default", exact = TRUE)))
    attrs$default <- .conductance_model_to_external(attr(conductance_model, "default", exact = TRUE),
                                                    conductance_model)
  if (!is.null(attr(conductance_model, "lower", exact = TRUE)))
    attrs$lower <- .conductance_model_to_external(attr(conductance_model, "lower", exact = TRUE),
                                                  conductance_model)
  if (!is.null(attr(conductance_model, "upper", exact = TRUE)))
    attrs$upper <- .conductance_model_to_external(attr(conductance_model, "upper", exact = TRUE),
                                                  conductance_model)
  for (nm in names(attrs))
    attr(external_model, nm) <- attrs[[nm]]

  info <- attr(external_model, "gaussian_scale_info", exact = TRUE)
  if (!is.null(info))
    attr(external_model, "gaussian_scale_info") <- info

  external_model
}

.run_terradish_optimizer <- function(theta, optfn, optimizer, control,
                                     lower = rep(-Inf, length(theta)),
                                     upper = rep(Inf, length(theta)))
{
  if (optimizer == "newton")
    return(BoxConstrainedNewton(theta, optfn,
                                lower = lower,
                                upper = upper,
                                control = control))
  BoxConstrainedBFGS(theta, optfn,
                     lower = lower,
                     upper = upper,
                     control = control)
}

.conductance_model_factory_attr <- function(model, attr_name, default = NULL)
{
  out <- attr(model, attr_name, exact = TRUE)
  if (is.null(out))
    default
  else
    out
}

.resolve_terradish_optimizer <- function(optimizer, n_theta,
                                         conductance_model_factory = NULL)
{
  preferred <- .conductance_model_factory_attr(conductance_model_factory,
                                               "preferred_optimizer")
  if (identical(optimizer, "auto") && !is.null(preferred))
    return(preferred)
  if (identical(optimizer, "auto") && as.integer(n_theta) > 3L)
    return("bfgs")
  if (identical(optimizer, "auto"))
    return("newton")
  optimizer
}

#' Optimize a parameterized conductance surface
#'
#' Uses maximum likelihood to fit a parameterized conductance surface to genetic data,
#' given a function relating spatial data to conductance (a "conductance model")
#' and a function relating resistance distance (covariance) to genetic distance
#' (a "measurement model").
#'
#' @param formula A formula with a matrix of observed genetic distances on the lhs, and covariates used in the creation of \code{data} on the rhs
#' @param data An object of class \code{terradish_graph} (see
#'   \code{\link{conductance_surface}})
#' @param conductance_model A function of class
#'   \code{terradish_conductance_model_factory} (see
#'   \code{\link{terradish_conductance_model_factory}})
#' @param measurement_model A function of class
#'   \code{terradish_measurement_model} (see
#'   \code{\link{terradish_measurement_model}})
#' @param nu Effective Wishart degrees of freedom passed to measurement models
#'   that require it, such as \code{\link{generalized_wishart}} and
#'   \code{\link{wishart_covariance}}.  For biallelic SNPs this is usually the
#'   number of retained SNPs.  For microsatellites, use the independent
#'   allele-frequency count, approximately \eqn{\sum_l (K_l - 1)} where
#'   \eqn{K_l} is the number of observed alleles at locus \eqn{l}.  Ignored by
#'   non-Wishart measurement models.
#' @param theta Starting values for optimization
#' @param leverage Compute influence measures and leverage?
#' @param nonnegative Force regression-like \code{measurement_model} to have nonnegative slope?
#' @param conductance Retained for backward compatibility. Only
#'   \code{conductance = TRUE} is currently implemented.
#' @param optimizer The optimization algorithm to use: \code{newton} uses the
#'   exact Hessian during optimization, with computational cost that grows
#'   linearly with the number of parameters; \code{bfgs} uses an approximate
#'   Hessian during optimization with cheaper iterations but often more steps;
#'   and \code{auto} selects \code{bfgs} when there are more than three
#'   conductance parameters and \code{newton} otherwise. Regardless of the
#'   optimizer, \code{terradish()} evaluates the fitted model once more at the
#'   final parameter values with the exact Hessian so standard errors and
#'   summaries use the same final derivative machinery.
#' @param control A list containing options for the optimization routine (see \code{\link{NewtonRaphsonControl}} for list)
#' @param validate Numerical validation of leverage via package \code{numDeriv} (very slow, use for debugging small examples)
#' @param cores Number of worker processes to use for Hessian and leverage calculations. \code{1} evaluates serially.
#' @param curvature Curvature used for optimization steps and for the returned
#'   covariance matrix. \code{"exact"} (default) uses the exact Hessian.
#'   \code{"gauss_newton"} uses the Gauss-Newton / Fisher-information
#'   approximation, which equals the Fisher information at the optimum (where it
#'   is positive semidefinite, giving a well-defined \code{vcov}), requires only
#'   first derivatives of the conductance model, and equals the exact Hessian at
#'   a well-fitting optimum.
#'   Standard errors from \code{summary()} are then the asymptotic
#'   information-based errors. With \code{leverage = TRUE} the leverage
#'   diagnostics inherit the same approximation.
#' @param solver Linear-system solver used for the reduced Laplacian. \code{"direct"} uses sparse Cholesky updates; \code{"auto"} conservatively chooses between the direct and AMG backends based on graph size and right-hand-side count; \code{"amg"} uses smoothed-aggregation AMG-preconditioned conjugate gradients; \code{"pcg"} uses incomplete-Cholesky preconditioned conjugate gradients; \code{"pcg_jacobi"} keeps the older Jacobi-preconditioned prototype.
#' @param solver_control Optional named list of solver settings passed to
#'   \code{\link{terradish_algorithm}}. For \code{solver = "direct"}, supported
#'   entries include \code{factorization}, \code{supernodal_min_vertices},
#'   \code{supernodal_max_rhs}, \code{perm}, and \code{solve_backend}; the
#'   latter can be \code{"matrix"} or the experimental \code{"cholmod_cpp"} and
#'   \code{"cholmod_cpp_cached"} backends.
#'   For \code{solver = "auto"}, supported selection entries include
#'   \code{auto_direct_max_vertices},
#'   \code{auto_amg_min_vertices}, and \code{auto_direct_max_rhs}. For
#'   \code{solver = "amg"} or \code{"auto"}, \code{terradish()} also
#'   understands an adaptive schedule with entries such as \code{adaptive},
#'   \code{tol_early}, \code{tol_mid}, \code{tol_final}, \code{maxit_early},
#'   \code{maxit_mid}, \code{maxit_final}, and \code{warmup_evals}. AMG controls
#'   also support \code{reuse_preconditioner} to reuse multigrid hierarchy
#'   information across nearby optimization steps when possible, and
#'   \code{reuse_preconditioner_max_age} to periodically rebuild a reused
#'   hierarchy. Direct supernodal factorizations can benefit from a threaded
#'   BLAS, but those thread counts are controlled by the external R/BLAS build
#'   rather than by \code{terradish()} itself.
#'   If a fit stops with a CHOLMOD or non-positive-definite reduced Laplacian
#'   message, the current conductance values made the graph numerically
#'   singular. Common fixes are to scale raster covariates, start conductance
#'   parameters closer to zero, check retained raster cells for missing or
#'   infinite values, and retry with \code{solver = "auto"} or
#'   \code{solver = "amg"} while diagnosing the problem.
#' @param approximation Exploratory approximation used during optimization.
#'   \code{"none"} uses the full focal set throughout. \code{"landmark"}
#'   optimizes first on a space-filling subset of focal populations and then
#'   refines on the full likelihood when supported by the measurement model.
#'   \code{"coarse_raster"} optimizes first on an aggregated raster and then
#'   optionally refines on the full-resolution graph. If multiple coarse
#'   factors are supplied, they are evaluated from coarsest to finest before the
#'   final full-resolution stage. This is an opt-in warm-start strategy, not a
#'   replacement for the exact full-resolution likelihood unless
#'   \code{exact_refine = FALSE}.
#' @param approximation_control Optional named list controlling the landmark or
#'   coarse-raster approximation. Landmark entries include
#'   \code{n_landmarks}, \code{fraction}, \code{min_landmarks},
#'   \code{max_landmarks}, \code{method} (\code{"spacefill"},
#'   \code{"random"}, or \code{"sequential"}), \code{seed}, and
#'   \code{exact_refine}. For \code{approximation = "coarse_raster"},
#'   supported entries include \code{factor}, \code{aggregate_fun},
#'   \code{directions}, \code{exact_refine}, and \code{refine_control}; this
#'   path requires \code{conductance_surface(..., saveStack = TRUE)}.
#'   \code{factor} may be a
#'   single positive integer or a vector such as \code{c(8, 4, 2)} for a
#'   multilevel coarse-to-fine warm start. The final reported fit is still
#'   evaluated on the full data. \code{refine_control} may be a
#'   \code{\link{NewtonRaphsonControl}} list used only for that final exact
#'   full-resolution refinement; by default the main \code{control} is reused.
#'   For
#'   \code{measurement_model = leastsquares}, landmark exact-refinement is
#'   guarded and the approximation stage is skipped.
#'
#' @details 
#' \figure{terradish-sticker.png}{options: style='float: right; width: 150px; margin-left: 12px;'}
#' By "parameterized conductance surface", what is meant is a model
#' where the per-vertex conductance (and thus resistance distance) is a function of
#' spatial covariates. The choice of function is referred to in this package as
#' the "conductance model". The inverse problem (and the purpose of this
#' package) is to estimate the parameters of the conductance model, by relating
#' the (unknown, modeled) resistance distance to observed genetic dissimilarity
#' via a probability model (referred to as the "measurement model" throughout
#' this package).
#'
#'
#' For example, a log-linear choice of conductance model is:
#'
#'   \code{vertex_conductance[i] = exp(covariates[i,] \%*\% theta)}
#'
#' where \code{theta} are unknown parameters, \code{covariates} is a design
#' matrix with each row representing a vertex in the graph (cell in the
#' raster). The conductance of a given edge -- an offdiagonal entry in the
#' graph Laplacian -- is:
#'
#'   \code{edge_conductance[i,j] = vertex_conductance[i] + vertex_conductance[j]}
#'
#' if argument \code{conductance = TRUE}, and
#'
#'   \code{edge_conductance[i,j] = 1/(1/vertex_conductance[i] + 1/vertex_conductance[j])}
#'
#' otherwise. This alternative parameterization is not currently implemented,
#' so \code{conductance} must be \code{TRUE}. This differs by a factor of two
#' from CIRCUITSCAPE, where the edge 
#' conductance/resistance is the average of the vertex conductance/resistance.
#'
#' \code{terradish} estimates \code{theta} (and thus the conductance) by maximum
#' likelihood; by finding the
#' values of \code{theta} (and associated conductance) that result in
#' resistance distances that are closest to the observed genetic distances,
#' according to some measure of fit (like least squares). The optimization is
#' done via Newton's method (default; requires computation of Hessian during
#' optimization), via the BFGS algorithm (requires gradient only during
#' optimization) if \code{optimizer = "bfgs"}, or via a simple parameter-count
#' heuristic if \code{optimizer = "auto"}. In all cases, the final returned fit
#' is evaluated with the exact Hessian at the optimized parameters for inference.
#'
#' For an explanation of how categorical spatial covariates are handled, see
#' \code{details} of \code{\link{conductance_surface}} and the examples below.
#' The dummy coding of categorical covariates is done by the function passed to
#' \code{conductance_model} (e.g.  \code{terradish::loglinear_conductance}).
#'
#' \strong{Large rasters.} For large continuous rasters, two opt-in helpers can
#' reduce runtime. First, build the graph with
#' \code{conductance_surface(..., crop_buffer = )} when focal sites occupy only
#' part of the raster. This removes vertices outside the buffered sampling
#' extent before fitting. Second, use
#' \code{approximation = "coarse_raster"} with
#' \code{approximation_control = list(factor = c(4, 2), exact_refine = TRUE)}
#' to optimize on one or more aggregated rasters before refining on the original
#' graph. With \code{exact_refine = TRUE}, the returned coefficients and
#' likelihood are from the full-resolution graph; the coarse fits are only
#' starting values. To keep the full-resolution cleanup deliberately short, pass
#' \code{refine_control = NewtonRaphsonControl(maxit = 2, ...)} inside
#' \code{approximation_control}. With \code{exact_refine = FALSE}, the result is
#' faster but approximate and should be interpreted as a screening fit.
#'
#' \strong{Gaussian scale-aware conductance.}
#' \code{\link{gaussian_smoothed_loglinear_conductance}} can also use
#' \code{approximation = "coarse_raster"}. In that case the Gaussian smoother is
#' rebuilt on each aggregated raster, while the full-resolution sigma bounds and
#' unit conversion are preserved so a warm-started \code{sigma} continues to
#' mean the same map-unit distance at every stage.
#'
#' Currently, all of the built-in conductance models in \code{terradish} use the
#' default contrast coding (e.g. for a categorical covariate with \code{K}
#' factors, the estimated parameters are the \code{K-1} mean differences
#' against a reference category). The intercept is excluded if it is not
#' identifiable.
#'
#' If the fit is on the boundary (e.g. no spatial genetic structure) or is the
#' null model of isolation-by-distance, the fitted object will not contain
#' influence/leverage/gradient/hessian.
#'
#' @return An object of class \code{terradish} containing the fitted conductance
#'   parameters, optimized nuisance parameters, log-likelihood, model
#'   comparison statistics, timing/evaluation diagnostics, and optional
#'   leverage diagnostics. See \code{\link{terradish_methods}} for the
#'   available S3 methods.
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
#' # scaling spatial covariates helps avoid numeric overflow
#' covariates <- c(terra::scale(melip.altitude),
#'                 terra::scale(melip.forestcover))
#' names(covariates) <- c("altitude", "forestcover")
#' 
#' # create parameterized conductance surface
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#' 
#' fit_nnls <- terradish(melip.Fst ~ altitude * forestcover, data = surface, 
#'                    terradish::loglinear_conductance, terradish::leastsquares)
#' summary(fit_nnls)
#' 
#' # a different "measurement_model" that incorporates dependence
#' # among pairwise measurements
#' fit_mlpe <- terradish(melip.Fst ~ altitude * forestcover, data = surface, 
#'                    terradish::loglinear_conductance, terradish::mlpe)
#' summary(fit_mlpe)
#'
#' # conductance surface with 95% CI
#' fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)
#' 
#' # test for an interaction using a likelihood ratio test
#' fit_mlpe_interaction <- terradish(melip.Fst ~ forestcover * altitude, data = surface, 
#'                                terradish::loglinear_conductance, terradish::mlpe)
#' anova(fit_mlpe, fit_mlpe_interaction)
#' 
#' # test against null model of IBD using a LRT
#' fit_mlpe_ibd <- terradish(melip.Fst ~ 1, data = surface, 
#'                        terradish::loglinear_conductance, terradish::mlpe)
#' anova(fit_mlpe, fit_mlpe_ibd)
#'
#' # categorical covariates:
#' # categorical raster layers should be factor-valued, see ?terra::as.factor
#' # and 'details' section of ?conductance_surface
#' forestcover_class <- cut(terra::values(melip.forestcover)[,1], breaks = c(0, 1/6, 1/3, 1))
#' melip.forestcover_cat <- terra::setValues(melip.forestcover, as.numeric(forestcover_class))
#' melip.forestcover_cat <- terra::as.factor(melip.forestcover_cat)
#'
#' RAT <- levels(melip.forestcover_cat)[[1]]
#' RAT$VALUE <- levels(forestcover_class) #explicitly define level names
#' levels(melip.forestcover_cat) <- RAT
#' 
#' covariates_cat <- c(melip.forestcover_cat, melip.altitude)
#' names(covariates_cat) <- c("forestcover", "altitude")
#' 
#' surface_cat <- conductance_surface(covariates_cat, melip.coords, directions = 8)
#' 
#' fit_mlpe_cat <- terradish(melip.Fst ~ forestcover + altitude, surface_cat, 
#'                        terradish::loglinear_conductance, terradish::mlpe)
#' summary(fit_mlpe_cat)
#'
#' @export

terradish <- function(formula, 
                   data,
                   conductance_model = loglinear_conductance, 
                   measurement_model = mlpe, 
                   nu = NULL, 
                   theta = NULL,
                   leverage = FALSE, 
                   nonnegative = TRUE, 
                   conductance = TRUE, 
                   optimizer = c("newton", "bfgs", "auto"), 
                   control = NewtonRaphsonControl(verbose = TRUE, ctol = 1e-6, ftol = 1e-6), 
                   validate = FALSE,
                   cores = 1L,
                   curvature = c("exact", "gauss_newton"),
                   solver = c("direct", "auto", "amg", "pcg", "pcg_jacobi", "block_cg"),
                   solver_control = NULL,
                   approximation = c("none", "landmark", "coarse_raster"),
                   approximation_control = NULL)
{
  stopifnot(inherits(formula, "formula"))
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  stopifnot(inherits(conductance_model, c("terradish_conductance_model_factory",
                                          "radish_conductance_model_factory")))
  stopifnot(inherits(measurement_model, c("terradish_measurement_model",
                                          "radish_measurement_model")))
  stopifnot(length(cores) == 1, is.numeric(cores), cores >= 1)
  if (!isTRUE(conductance))
    stop("`conductance = FALSE` is not currently supported.", call. = FALSE)
  solver <- match.arg(solver)
  curvature <- match.arg(curvature)
  approximation <- match.arg(approximation)

  # get response, remove lhs from formula
  terms    <- terms(formula)
  vars     <- as.character(attr(terms, "variables"))[-1]
  response <- attr(terms, "response")
  S        <- if(response) eval(attr(terms, "variables")[[response + 1L]], parent.frame())
              else stop("'formula' must have genetic distance matrix on lhs")
  is_ibd   <- length(vars) == 1
  formula  <- if (!is_ibd) reformulate(attr(terms, "term.labels"))
              else formula(~1)

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model_factory <- conductance_model
  rebuild_conductance_model_for_surface <- .conductance_model_factory_attr(
    conductance_model_factory,
    "rebuild_for_surface",
    NULL
  )
  if (identical(approximation, "coarse_raster") &&
      isTRUE(.conductance_model_factory_attr(conductance_model_factory,
                                             "requires_fixed_graph", FALSE)) &&
      is.null(rebuild_conductance_model_for_surface))
    stop("`approximation = \"coarse_raster\"` is not currently supported for this conductance model")
  conductance_model <- conductance_model_factory(formula, data$x)
  conductance_model_user <- .externalize_conductance_model(conductance_model)
  conductance_supports_partial <- !identical(
    attr(conductance_model, "supports_partial", exact = TRUE),
    FALSE
  )

  # initialize theta
  default <- attr(conductance_model, "default")
  bounds <- .conductance_model_bounds(conductance_model)
  if (is.null(theta))
    theta <- default
  else
  {
    stopifnot(length(theta) == length(default))
    names(theta) <- names(.conductance_model_parameter_scale(conductance_model))
    theta <- .conductance_model_to_internal(theta, conductance_model)
  }
  if (any(theta < bounds$lower | theta > bounds$upper))
    stop("Starting values in `theta` must lie within the conductance-model bounds",
         call. = FALSE)

  optimizer <- .resolve_terradish_optimizer(match.arg(optimizer), length(theta),
                                           conductance_model_factory = conductance_model_factory)
  if (isTRUE(leverage) && !conductance_supports_partial)
  {
    if (isTRUE(attr(conductance_model, "gaussian_scale", exact = TRUE)))
      warning("Leverage diagnostics are not available for Gaussian scale-aware conductance because `partial_X` currently assumes a one-to-one mapping between conductance parameters and design-matrix covariates; models with separate sigma parameters do not satisfy that contract.",
              call. = FALSE)
    else
      warning("Leverage diagnostics are not available for this conductance model; skipping `partial_X` calculations.",
              call. = FALSE)
    leverage <- FALSE
  }
  fcalls    <- new("FunctionCall", count = 0L)
  diagnostics <- .terradish_new_diagnostics()
  control$diagnostics <- diagnostics
  make_optfn <- function(eval_data,
                         eval_S,
                         phi_state,
                         solver_state,
                         eval_conductance_model = conductance_model)
  {
    function(par, gradient, hessian)
    {
      fcalls$count <- fcalls$count + 1L
      current_solver_control <- .terradish_solver_control_for_phase(
        solver = solver,
        solver_control = solver_control,
        eval_count = fcalls$count,
        final = FALSE
      )
      fit <- terradish_algorithm(f = eval_conductance_model,
                              g = measurement_model,
                              s = eval_data,
                              S = eval_S,
                              nu = nu,
                              phi = phi_state$value,
                              theta = c(par),
                              gradient = gradient,
                              hessian = hessian,
                              partial = FALSE,
                              nonnegative = nonnegative,
                              cores = cores,
                              curvature = curvature,
                              solver = solver,
                              solver_control = current_solver_control,
                              solver_warm_start = solver_state$warm_start,
                              solver_reuse_state = solver_state$reuse_state)
      phi_state$value <- fit$phi
      solver_state$warm_start <- fit$solver_warm_start
      solver_state$reuse_state <- fit$solver_reuse_state
      .terradish_record_algorithm_diagnostics(
        diagnostics,
        fit,
        objective = TRUE,
        gradient = gradient,
        hessian = hessian,
        partial = FALSE
      )
      fit
    }
  }

  approximation_info <- list(type = "none", used = FALSE)

  if (!is_ibd)
  {
    measurement_model_name <- .terradish_measurement_model_name(measurement_model)
    landmark_problem <- NULL
    coarse_problem <- NULL
    if (identical(approximation, "landmark"))
    {
      landmark_problem <- .terradish_landmark_subset(data, S, approximation_control = approximation_control)
      if (isTRUE(landmark_problem$used) &&
          isTRUE(landmark_problem$control$exact_refine) &&
          identical(measurement_model_name, "leastsquares"))
      {
        approximation_info <- list(
          type = "landmark",
          used = FALSE,
          stage = "disabled_for_leastsquares_exact_refine",
          n_landmarks = landmark_problem$control$n_landmarks,
          full_focal = nrow(S),
          focal_fraction = landmark_problem$control$n_landmarks / nrow(S),
          method = landmark_problem$control$method,
          landmark_index = landmark_problem$index,
          exact_refine = TRUE,
          refine_guard = "disabled_for_leastsquares"
        )
        landmark_problem <- NULL
      }
    }
    else if (identical(approximation, "coarse_raster"))
    {
      coarse_problem <- .coarse_raster_surface(
        data,
        approximation_control = approximation_control
      )
      approximation_info <- list(
        type = "coarse_raster",
        used = isTRUE(coarse_problem$used),
        stage = if (isTRUE(coarse_problem$used) &&
                    isTRUE(coarse_problem$control$exact_refine))
          "coarse_then_exact_refine"
        else if (isTRUE(coarse_problem$used))
          "optimization_only"
        else
          "skipped_factor_1",
        factor = coarse_problem$control$factor,
        n_levels = length(coarse_problem$stages),
        directions = coarse_problem$control$directions,
        full_vertices = coarse_problem$full_vertices,
        coarse_vertices = coarse_problem$coarse_vertices,
        full_focal = nrow(S),
        duplicate_demes = coarse_problem$duplicate_demes,
        unique_demes = coarse_problem$unique_demes,
        exact_refine = isTRUE(coarse_problem$control$exact_refine)
      )
    }

    iters <- 0L
    exact_phi_state <- new.env(parent = emptyenv())
    exact_phi_state$value <- NULL
    exact_solver_state <- new.env(parent = emptyenv())
    exact_solver_state$warm_start <- NULL
    exact_solver_state$reuse_state <- NULL
    exact_control <- control
    exact_control_source <- "shared"
    coarse_stage_details <- NULL
    exact_refine_steps <- NA_integer_

    if (!is.null(landmark_problem) && isTRUE(landmark_problem$used))
    {
      approximation_info <- list(
        type = "landmark",
        used = TRUE,
        stage = "optimization_only",
        n_landmarks = landmark_problem$control$n_landmarks,
        full_focal = nrow(S),
        focal_fraction = landmark_problem$control$n_landmarks / nrow(S),
        method = landmark_problem$control$method,
        landmark_index = landmark_problem$index,
        exact_refine = isTRUE(landmark_problem$control$exact_refine)
      )

      approx_phi_state <- new.env(parent = emptyenv())
      approx_phi_state$value <- NULL
      approx_solver_state <- new.env(parent = emptyenv())
      approx_solver_state$warm_start <- NULL
      approx_solver_state$reuse_state <- NULL

      approx_problem <- .run_terradish_optimizer(
        theta = theta,
        optfn = make_optfn(landmark_problem$data, landmark_problem$S, approx_phi_state, approx_solver_state),
        optimizer = optimizer,
        control = control,
        lower = bounds$lower,
        upper = bounds$upper
      )
      iters <- iters + approx_problem$iters
      theta <- c(approx_problem$par)
      names(theta) <- names(default)

      exact_phi_state$value <- approx_phi_state$value
      exact_solver_state$reuse_state <- approx_solver_state$reuse_state
    }
    else if (!is.null(coarse_problem) && isTRUE(coarse_problem$used))
    {
      approx_phi_state <- new.env(parent = emptyenv())
      approx_phi_state$value <- NULL
      coarse_stage_details <- vector("list", length(coarse_problem$stages))
      names(coarse_stage_details) <- names(coarse_problem$stages)

      for (stage_name in names(coarse_problem$stages))
      {
        stage <- coarse_problem$stages[[stage_name]]
        coarse_conductance_model <- if (!is.null(rebuild_conductance_model_for_surface))
          rebuild_conductance_model_for_surface(
            formula,
            stage$surface,
            reference_model = conductance_model
          )
        else
          conductance_model_factory(formula, stage$surface$x)
        coarse_bounds <- .conductance_model_bounds(coarse_conductance_model)
        approx_solver_state <- new.env(parent = emptyenv())
        approx_solver_state$warm_start <- NULL
        approx_solver_state$reuse_state <- NULL

        approx_problem <- .run_terradish_optimizer(
          theta = theta,
          optfn = make_optfn(stage$surface,
                             S,
                             approx_phi_state,
                             approx_solver_state,
                             eval_conductance_model = coarse_conductance_model),
          optimizer = optimizer,
          control = control,
          lower = coarse_bounds$lower,
          upper = coarse_bounds$upper
        )
        iters <- iters + approx_problem$iters
        theta <- c(approx_problem$par)
        names(theta) <- names(default)
        coarse_stage_details[[stage_name]] <- list(
          factor = stage$factor,
          vertices = stage$coarse_vertices,
          duplicate_demes = stage$duplicate_demes,
          unique_demes = stage$unique_demes,
          optimizer_steps = approx_problem$iters
        )
      }

      # The coarse graph has a different dimension, so only transfer the
      # nuisance-parameter warm start to the exact full-resolution stage.
      exact_phi_state$value <- approx_phi_state$value
      if (isTRUE(coarse_problem$control$exact_refine) &&
          !is.null(coarse_problem$control$refine_control))
      {
        exact_control <- coarse_problem$control$refine_control
        exact_control$diagnostics <- diagnostics
        exact_control_source <- "custom"
      }
    }

    if (!isTRUE(approximation_info$used) || isTRUE(approximation_info$exact_refine))
    {
      exact_problem <- .run_terradish_optimizer(
        theta = theta,
        optfn = make_optfn(data, S, exact_phi_state, exact_solver_state),
        optimizer = optimizer,
        control = exact_control,
        lower = bounds$lower,
        upper = bounds$upper
      )
      iters <- iters + exact_problem$iters
      theta <- c(exact_problem$par)
      names(theta) <- names(default)
      if (!is.null(coarse_problem) && isTRUE(coarse_problem$used))
        exact_refine_steps <- exact_problem$iters
    }

    if (identical(approximation_info$type, "coarse_raster") &&
        isTRUE(approximation_info$used))
    {
      approximation_info$stages <- coarse_stage_details
      approximation_info$coarse_steps <- sum(vapply(coarse_stage_details,
                                                    `[[`,
                                                    integer(1),
                                                    "optimizer_steps"))
      approximation_info$refine_steps <- if (isTRUE(approximation_info$exact_refine))
        exact_refine_steps
      else
        0L
      approximation_info$refine_control <- if (isTRUE(approximation_info$exact_refine))
        exact_control_source
      else
        "not_used"
    }
  }
  else
  {
    iters <- 0L
    theta <- c(default)
    exact_phi_state <- new.env(parent = emptyenv())
    exact_phi_state$value <- NULL
    exact_solver_state <- new.env(parent = emptyenv())
    exact_solver_state$warm_start <- NULL
    exact_solver_state$reuse_state <- NULL
  }

  final_solver_control <- .terradish_solver_control_for_phase(
    solver = solver,
    solver_control = solver_control,
    eval_count = fcalls$count,
    final = TRUE
  )
  fit <- terradish_algorithm(f = conductance_model, g = measurement_model, 
                          s = data, S = S, nu = nu, theta = theta, phi = exact_phi_state$value,
                          gradient = TRUE, hessian = TRUE, partial = leverage,
                          nonnegative = nonnegative, cores = cores,
                          curvature = curvature,
                          solver = solver, solver_control = final_solver_control,
                          solver_warm_start = exact_solver_state$warm_start,
                          solver_reuse_state = exact_solver_state$reuse_state)
  .terradish_record_algorithm_diagnostics(
    diagnostics,
    fit,
    objective = TRUE,
    gradient = TRUE,
    hessian = TRUE,
    partial = leverage
  )
  fit$solver_warm_start <- NULL
  fit$solver_reuse_state <- NULL

  fit$response <- S
  fit$gradient_internal <- fit$gradient
  fit$hessian_internal <- fit$hessian
  fit$gradient <- .conductance_model_gradient_to_external(fit$gradient, conductance_model)
  fit$hessian <- .conductance_model_hessian_to_external(fit$hessian, conductance_model)
  theta_external <- .conductance_model_to_external(theta, conductance_model)
  # attach parameter names so users can index hessian by name
  if (!is.null(names(theta_external)))
    dimnames(fit$hessian) <- list(names(theta_external), names(theta_external))

  if (fit$boundary)
    warning("Optimum for subproblem is on boundary (e.g. no spatial genetic structure): cannot optimize theta.\nTry different starting values.")
  no_coef <- fit$boundary || is_ibd 

  # calculate leverage for genetic distance and spatial covariates
  leverage <- leverage && !no_coef
  if (leverage)
  {
    ihess      <- .safe_hessian_inverse(fit$hessian)
    leverage_S <- -matrix(fit$partial_S, length(S), length(theta_external)) %*% ihess
    leverage_S <- array(leverage_S, dim = c(nrow(S), ncol(S), length(theta_external)))
    for (k in seq_along(theta_external))
      leverage_S[, , k] <- (leverage_S[, , k] + t(leverage_S[, , k])) / 2
    leverage_X <- array(NA, dim = dim(fit$partial_X))
    for (k in 1:length(theta_external))
      leverage_X[,,k] <- -fit$partial_X[,,k] %*% ihess
  }
  num_leverage_S <- NULL
  num_leverage_X <- NULL

  # numerical validation of derivatives
  validate <- validate && !no_coef
  if (validate)
  {
    warning("`validate = TRUE` is not currently implemented for `terradish()`; skipping numerical leverage validation.",
            call. = FALSE)
    validate <- FALSE
  }

  out <- list(call           = match.call(),
              formula        = formula,
              dim            = c("vertices" = nrow(data$x),
                                 "focal"    = nrow(S),
                                 "edge"     = ncol(data$adj)),
              cost           = c("newton_steps"   = iters,
                                 "function_calls" = fcalls$count + 1),
              diagnostics    = .terradish_diagnostics_snapshot(diagnostics),
              submodels      = list("f" = conductance_model_user,
                                    "f_internal" = conductance_model,
                                    "f_factory" = conductance_model_factory,
                                    "g" = measurement_model),
              fit            = fit,
              loglik         = -fit$objective,
              df             = (!no_coef)*length(theta_external) + length(fit$phi),
              aic            = 2*fit$objective + 2*(!no_coef)*length(theta_external) + 2*length(fit$phi),
              mle            = list("theta"    = if(no_coef) NULL else theta_external,
                                    "theta_internal" = if(no_coef) NULL else theta,
                                    "gradient" = if(no_coef) NULL else -fit$gradient,
                                    "gradient_internal" = if(no_coef) NULL else -fit$gradient_internal,
                                    "hessian"  = if(no_coef) NULL else -fit$hessian,
                                    "hessian_internal" = if(no_coef) NULL else -fit$hessian_internal),
              approximation  = approximation_info,
              leverage       = list("S" = if(!leverage) NULL else leverage_S,
                                    "X" = if(!leverage) NULL else leverage_X,
                                    "validate" = if(!validate || !leverage) NULL 
                                                 else list("S" = num_leverage_S,
                                                           "X" = num_leverage_X))
              )
  class(out) <- c("terradish", "radish")
  out
}

#' Legacy radish fit wrapper
#'
#' Deprecated compatibility wrapper retained for older code that still calls
#' \code{radish()}.
#'
#' @param ... Arguments passed through the deprecated \code{radish()}
#'   compatibility wrapper to \code{\link{terradish}}.
#' @name legacy_radish_fit_wrapper
#' @keywords internal
NULL

#' @rdname legacy_radish_fit_wrapper
#' @export
radish <- function(...)
{
  .terradish_deprecate("radish", "terradish")
  .terradish_forward_call(match.call(), "terradish")
}

#' Methods for fitted terradish models
#'
#' S3 methods for working with fitted objects returned by \code{\link{terradish}}
#' and summaries returned by \code{\link[base:summary]{summary()}}.
#'
#' @param x A fitted \code{terradish} object.
#' @param object,alternative Fitted \code{terradish} objects.
#' @param digits Number of digits to print.
#' @param signif.stars Should significance stars be printed in coefficient
#'   tables?
#' @param type Type of fitted values to return.
#' @param nsim Number of simulated response matrices to generate.
#' @param seed Optional random seed used by \code{simulate()}.
#' @param k Penalty multiplier supplied to \code{AIC()}.
#' @param method Simulation method.
#' @param conf.level Confidence level for nuisance-parameter confidence
#'   intervals reported by \code{summary()}.  These intervals are conditional on
#'   the fitted conductance surface.
#' @param ... Additional arguments passed through to generic methods.
#'
#' @return
#' \itemize{
#'   \item \code{print()} returns its input invisibly.
#'   \item \code{summary()} returns an object of class
#'     \code{summary.terradish}, with legacy \code{summary.radish}
#'     compatibility retained.
#'   \item \code{coef()} returns the fitted conductance coefficients.
#'   \item \code{fitted()} returns fitted responses, distances, or covariance.
#'   \item \code{simulate()} returns one or more simulated response matrices.
#'   \item \code{anova()} returns a likelihood-ratio comparison table.
#'   \item \code{logLik()} returns a \code{logLik} object.
#'   \item \code{AIC()} returns Akaike's Information Criterion.
#'   \item \code{residuals()} returns the residual genetic-distance matrix.
#' }
#'
#' @name terradish_methods
NULL

#' Legacy radish S3 methods
#'
#' Deprecated S3 methods retained for legacy objects that still dispatch on the
#' older \code{"radish"} class names.
#'
#' @param x A fitted legacy \code{radish} object.
#' @param object,alternative Fitted legacy \code{radish} objects.
#' @param digits Number of digits to print.
#' @param signif.stars Should significance stars be printed in coefficient
#'   tables?
#' @param type Type of fitted values to return.
#' @param nsim Number of simulated response matrices to generate.
#' @param seed Optional random seed used by \code{simulate()}.
#' @param k Penalty multiplier supplied to \code{AIC()}.
#' @param method Simulation method.
#' @param conf.level Confidence level for nuisance-parameter confidence
#'   intervals reported by \code{summary()}. These intervals are conditional on
#'   the fitted conductance surface.
#' @param ... Additional arguments passed through to generic methods.
#' @name legacy_radish_methods
#' @keywords internal
NULL

#' @rdname legacy_radish_methods
#' @method print radish
#' @export
print.radish <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("Conductance surface estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (!x$fit$boundary && !is.null(x$mle$theta))
  {
    cat("Coefficients:\n")
    print.default(format(x$mle$theta, digits = digits), print.gap = 2L, quote = FALSE)
  }
  else
  {
    cat("No coefficients\n")
  }
  cat("\n")
  cat("Loglikelihood:", x$loglik, paste0("(", x$df), "degrees freedom)   AIC:", x$aic, "\n")
  invisible(x)
}

#' @rdname legacy_radish_methods
#' @method summary radish
#' @export
summary.radish <- function(object, conf.level = 0.95, ...)
{
  x <- object
  tol <- sqrt(.Machine$double.eps) #for checking singularity

  no_coef <- x$fit$boundary || is.null(x$mle$theta)
  if (!no_coef)
  {
    ztable <- matrix(0, length(x$mle$theta), 4)
    colnames(ztable)      <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(ztable)      <- names(x$mle$theta)
    ztable[,"Estimate"]   <- x$mle$theta
    ztable[,"Std. Error"] <- sqrt(diag(solve(x$fit$hessian)))
    ztable[,"z value"]    <- ztable[,"Estimate"]/ztable[,"Std. Error"]
    ztable[,"Pr(>|z|)"]   <- pmin(2*(1 - pnorm(abs(ztable[,"z value"]))), 1)

    ehess <- eigen(x$fit$hessian)

    if (any(ehess$values < 0))
      warning("Hessian matrix has negative eigenvalues: possibly a saddle point")
    if (any(abs(ehess$values) < tol * max(abs(ehess$values))))
      warning("Hessian matrix is singular or nearly singular: model is probably non-identifiable")

    vcov <- ehess$vectors %*% diag(1/ehess$values, nrow = nrow(x$fit$hessian)) %*% t(ehess$vectors)
    vcor <- cov2cor(vcov)
    rownames(vcor) <- colnames(vcor) <- 
      rownames(vcov) <- colnames(vcov) <- 
        rownames(ztable)
  }

  out <- list(boundary      = x$fit$boundary,
              phi           = x$fit$phi[,1],
              phi_table     = NULL,
              phi_vcov      = NULL,
              ztable        = if (no_coef) NULL else ztable,
              vcor          = if (no_coef) NULL else vcor,
              vcov          = if (no_coef) NULL else vcov,
              gradnorm      = if (no_coef) NA else sqrt(sum(x$mle$gradient^2)),
              loglik        = x$loglik,
              df            = x$df,
              aic           = x$aic,
              fcalls        = x$cost["function_calls"],
              iters         = x$cost["newton_steps"],
              diagnostics   = x$diagnostics,
              call          = x$call,
              dim           = x$dim
              )

  phi_summary <- .conditional_phi_table(x$fit$phi[,1], x$fit$phi_hessian,
                                        conf.level = conf.level)
  if (!is.null(phi_summary))
  {
    out$phi_table <- phi_summary$table
    out$phi_vcov <- phi_summary$vcov
  }

  class(out) <- c("summary.terradish", "summary.radish")
  out
}

#' @rdname legacy_radish_methods
#' @method print summary.radish
#' @export
print.summary.radish <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...)
{
  cat("Conductance surface with", x$dim["vertices"], "vertices", 
      paste0("(", x$dim["focal"]), "focal) estimated by maximum likelihood\n")
  cat("Call:   ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Loglikelihood:", x$loglik, paste0("(", x$df), "degrees freedom)\nAIC:", x$aic, "\n\n")
  cat("Number of function calls:", x$fcalls, "\n")
  cat("Number of optimization steps:", x$iters, "\n")
  cat("Norm of gradient at MLE:", x$gradnorm, "\n\n")
  if (length(x$phi))
  {
    cat("Nuisance parameters")
    if (!is.null(x$phi_table))
      cat(" (conditional on fitted conductance surface)")
    cat(":\n")
    if (!is.null(x$phi_table))
      print.default(format(x$phi_table, digits = digits), print.gap = 2L, quote = FALSE)
    else
      print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  if (!x$boundary && !is.null(x$ztable))
  {
    cat("Coefficients:\n")
    printCoefmat(x$ztable, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    if (nrow(x$ztable) > 1)
    {
      cat("\n")
      cat("Correlation of Coefficients:\n")
      print(as.dist(x$vcor))
    }
  }
  else if (x$boundary)
  {
    cat("Model fit is on boundary (e.g. no genetic structure), no meaningful coefficients\n")
  }
  else
  {
    cat("No coefficients\n")
  }
  invisible(x)
}

#' @rdname legacy_radish_methods
#' @method coef radish
#' @export
coef.radish <- function(object, ...)
{
  if (!is.null(object$mle$theta))
    return(object$mle$theta)
  else
    return(c())
}

#' @rdname legacy_radish_methods
#' @method fitted radish
#' @export
fitted.radish <- function(object, type = c("response", "distance", "covariance"), ...)
{
  type <- match.arg(type)
  if (type == "response")
    as.matrix(object$fit$fitted)
  else if (type == "distance")
    dist_from_cov(as.matrix(object$fit$covariance))
  else if (type == "covariance")
    as.matrix(object$fit$covariance)
}

#' @rdname legacy_radish_methods
#' @method simulate radish
#' @export
simulate.radish <- function(object, nsim = 1, seed = NULL, method = c("permutation", "parametric"), ...)
{
  if (!is.null(seed))
    set.seed(seed)
  method <- match.arg(method)
  if (method == "parametric")
  {
    stop("Parametric simulation not yet supported")
  } 
  else if (method == "permutation") 
  {
    fit <- fitted(object)
    resid  <- residuals(object)
    sims   <- array(NA, c(nrow(resid), ncol(resid), nsim))
    for (i in 1:nsim)
    {
      ind <- sample(1:nrow(resid))
      sims[,,i] <- fit + resid[ind,ind]
    }
  }
  if (nsim == 1) sims[,,1] else sims
}

#' @rdname legacy_radish_methods
#' @method anova radish
#' @export
anova.radish <- function(object, ..., alternative = NULL)
{
  dots <- list(...)
  if (is.null(alternative))
    alternative <- dots[[1]]
  stopifnot(inherits(object, c("terradish", "radish")) &&
            inherits(alternative, c("terradish", "radish")))
  stopifnot(!object$fit$boundary && !alternative$fit$boundary)

  if (object$df >= alternative$df)
  {
    full <- object
    reduced <- alternative
  } 
  else
  {
    full <- alternative
    reduced <- object
  }

  form_reduced <- paste("Null:", paste(reduced$formula, collapse = " "))
  form_full    <- paste("Alt:", paste(full$formula, collapse = " "))

  Chisq <- 2 * (full$loglik - reduced$loglik)
  Df    <- full$df - reduced$df
  P     <- pchisq(Chisq, Df, lower.tail = FALSE)
  Ll    <- c(reduced$loglik, full$loglik)
  Np    <- c(reduced$df, full$df)

  out <- cbind("logLik" = Ll, 
               "Df" = Np, 
               "ChiSq" = c(NA, Chisq), 
               "Df(ChiSq)" = c(NA, Df), 
               "Pr(>Chi)" = c(NA, P))
  rownames(out) <- c("Null", "Alt")

  attr(out, "heading") <- c("Likelihood ratio test",
                           form_reduced, form_full)
  class(out) <- "anova"
  out
}

#' @rdname legacy_radish_methods
#' @method logLik radish
#' @export
logLik.radish <- function(object, ...)
{
  val <- object$loglik
  attr(val, "df") <- object$df
  class(val) <- "logLik"
  val
}

#' @rdname legacy_radish_methods
#' @method AIC radish
#' @export
AIC.radish <- function(object, ..., k = 2)
{
  object$aic
}

#' @rdname legacy_radish_methods
#' @method residuals radish
#' @export
residuals.radish <- function(object, ...)
{
  fit <- fitted(object)
  object$fit$response - fit
}

#' @rdname terradish_methods
#' @method print terradish
#' @export
print.terradish <- print.radish

#' @rdname terradish_methods
#' @method summary terradish
#' @export
summary.terradish <- summary.radish

#' @rdname terradish_methods
#' @method print summary.terradish
#' @export
print.summary.terradish <- print.summary.radish

#' @rdname terradish_methods
#' @method coef terradish
#' @export
coef.terradish <- coef.radish

#' @rdname terradish_methods
#' @method fitted terradish
#' @export
fitted.terradish <- fitted.radish

#' @rdname terradish_methods
#' @method simulate terradish
#' @export
simulate.terradish <- simulate.radish

#' @rdname terradish_methods
#' @method anova terradish
#' @export
anova.terradish <- anova.radish

#' @rdname terradish_methods
#' @method logLik terradish
#' @export
logLik.terradish <- logLik.radish

#' @rdname terradish_methods
#' @method AIC terradish
#' @export
AIC.terradish <- AIC.radish

#' @rdname terradish_methods
#' @method residuals terradish
#' @export
residuals.terradish <- residuals.radish

# Internal compatibility aliases retained for housekeeping-sized refactors.
.radish_amg_control_defaults <- .terradish_amg_control_defaults
.radish_solver_control_for_phase <- .terradish_solver_control_for_phase
.radish_landmark_control_defaults <- .terradish_landmark_control_defaults
.radish_measurement_model_name <- .terradish_measurement_model_name
.radish_landmark_subset <- .terradish_landmark_subset
.run_radish_optimizer <- .run_terradish_optimizer
