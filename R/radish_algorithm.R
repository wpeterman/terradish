.terradish_algorithm_derivative_chunk <- function(idx, state)
{
  if (!length(idx))
    return(list())

  chunk_state <- vector("list", length(idx))
  rhs_blocks <- vector("list", length(idx))

  for (m in seq_along(idx))
  {
    k <- idx[m]
    dgrad__ddl_dC <- state$C$df__dtheta(k)
    dgrad__ddl_dQnG <- laplacian_derivative_matrix_product(
      dgrad__ddl_dC,
      state$s$adj,
      state$G
    )
    dgrad__ddl_dE <- -state$tG %*% dgrad__ddl_dQnG
    dgrad__dE <- state$subproblem$jacobian_E(dgrad__ddl_dE)

    chunk_state[[m]] <- list(k = k, dgrad__ddl_dE = dgrad__ddl_dE)
    rhs_blocks[[m]] <- graph_rhs_matrix_product(
      state$s$demes,
      state$N,
      dgrad__dE
    )
    if (!isTRUE(state$gauss_newton))
      rhs_blocks[[m]] <- rhs_blocks[[m]] - 2 * dgrad__ddl_dQnG %*% state$dl_dE
  }

  rhs_block <- if (length(rhs_blocks) == 1L) rhs_blocks[[1]] else do.call(cbind, rhs_blocks)
  solve_result <- .terradish_solver_solve(state$solver_state, rhs_block)
  solved_block <- as.matrix(solve_result$solution)
  n_rhs <- ncol(state$Zn)

  out <- lapply(seq_along(idx), function(m)
  {
    k <- chunk_state[[m]]$k
    cols <- ((m - 1L) * n_rhs + 1L):(m * n_rhs)
    dgrad__dQnG <- t(solved_block[, cols, drop = FALSE])
    dgrad__dC <- backpropagate_laplacian_to_conductance(dgrad__dQnG, state$tG, state$s$adj)

    hess_row <- c(crossprod(state$df__dtheta_matrix, c(dgrad__dC)))
    if (!isTRUE(state$gauss_newton))
      for (l in seq_along(state$theta))
        hess_row[l] <- hess_row[l] +
          c(state$dl_dC) %*% state$C$d2f__dtheta_dtheta(k, l)

    partial_X_k <- NULL
    partial_S_k <- NULL
    if (state$partial)
    {
      partial_X_k <- matrix(0, state$N, length(state$theta))
      for (l in seq_along(state$theta))
        partial_X_k[, l] <- c(dgrad__dC) * state$C$df__dx(l) +
          c(state$dl_dC) * state$C$d2f__dtheta_dx(k, l)
      partial_S_k <- state$subproblem$jacobian_S(chunk_state[[m]]$dgrad__ddl_dE)
    }

    list(k = k,
         hess_row = hess_row,
         partial_X_k = partial_X_k,
         partial_S_k = partial_S_k)
  })
  attr(out, "solver_info") <- solve_result$info
  out
}

.graph_reduced_index <- function(s, n_vertices)
{
  idx <- s$reduced_index
  if (is.null(idx))
    idx <- seq_len(n_vertices - 1L)
  idx
}

.graph_edge_pairs <- function(s, n_vertices)
{
  edge_pairs <- s$edge_pairs
  if (!is.null(edge_pairs))
  {
    edge_pairs <- as.matrix(edge_pairs)
    if (nrow(edge_pairs) == 2L && ncol(edge_pairs) != 2L)
      edge_pairs <- t(edge_pairs)
    storage.mode(edge_pairs) <- "integer"
    return(edge_pairs)
  }

  adj <- s$adj
  if (is.null(adj))
    stop("`terradish_graph` is missing both `edge_pairs` and `adj`")
  edge_pairs <- t(adj) + 1L
  storage.mode(edge_pairs) <- "integer"
  edge_pairs
}

.graph_rhs <- function(s, n_vertices)
{
  rhs <- s$rhs
  if (!is.null(rhs))
    return(rhs)

  rhs <- matrix(-1 / n_vertices, nrow = n_vertices - 1L, ncol = length(s$demes))
  keep <- s$demes < n_vertices
  if (any(keep))
  {
    idx <- cbind(s$demes[keep], which(keep))
    rhs[idx] <- rhs[idx] + 1
  }
  rhs
}

.conductance_df_dtheta_matrix <- function(C, theta)
{
  if (!is.null(C$df__dtheta_matrix))
    return(as.matrix(C$df__dtheta_matrix))

  vapply(
    seq_along(theta),
    function(k) c(C$df__dtheta(k)),
    numeric(length(C$conductance))
  )
}

.graph_reduced_laplacian <- function(s, conductance)
{
  edge_pairs <- .graph_edge_pairs(s, length(conductance))
  if (!is.null(edge_pairs))
    return(forceSymmetric(assemble_reduced_laplacian(conductance, edge_pairs)))

  n_vertices <- length(conductance)
  reduced_index <- .graph_reduced_index(s, n_vertices)
  Q <- s$laplacian
  Q@x[] <- -conductance[s$adj[1,] + 1] - conductance[s$adj[2,] + 1]
  Qd <- Diagonal(n_vertices, x = -rowSums(Q))
  forceSymmetric((Q + Qd)[reduced_index, reduced_index, drop = FALSE])
}

.terradish_auto_solver_defaults <- function()
{
  list(
    auto_direct_max_vertices = 750000L,
    auto_amg_min_vertices = 1500000L,
    auto_direct_max_rhs = 64L
  )
}

.terradish_resolve_solver <- function(s, solver, n_vertices, solver_control = NULL)
{
  rhs <- .graph_rhs(s, n_vertices)
  n_rhs <- ncol(rhs)

  if (!identical(solver, "auto"))
    return(list(type = solver,
                requested_type = solver,
                solver_control = solver_control,
                reason = NULL,
                n_vertices = as.integer(n_vertices),
                n_rhs = as.integer(n_rhs)))

  control_list <- if (is.null(solver_control)) list() else as.list(solver_control)
  auto_defaults <- .terradish_auto_solver_defaults()
  auto_overrides <- control_list[intersect(names(control_list), names(auto_defaults))]
  auto_control <- modifyList(auto_defaults, auto_overrides)

  if (n_vertices <= as.integer(auto_control$auto_direct_max_vertices))
  {
    resolved <- "direct"
    reason <- "graph_not_large_enough_for_amg"
  }
  else if (n_vertices < as.integer(auto_control$auto_amg_min_vertices))
  {
    resolved <- "direct"
    reason <- "prefer_direct_until_larger_graphs"
  }
  else if (n_rhs <= as.integer(auto_control$auto_direct_max_rhs))
  {
    resolved <- "amg"
    reason <- "graph_large_enough_for_amg"
  }
  else
  {
    resolved <- "direct"
    reason <- "large_graph_with_many_rhs_favors_direct"
  }

  list(type = resolved,
       requested_type = "auto",
       solver_control = solver_control,
       reason = reason,
       n_vertices = as.integer(n_vertices),
       n_rhs = as.integer(n_rhs))
}

.terradish_amg_reuse_signature <- function(control)
{
  list(
    coarse_enough = as.integer(control$coarse_enough),
    npre = as.integer(control$npre),
    npost = as.integer(control$npost),
    sa_relax = control$sa_relax,
    aggr_eps_strong = control$aggr_eps_strong,
    estimate_spectral_radius = isTRUE(control$estimate_spectral_radius),
    power_iters = as.integer(control$power_iters)
  )
}

.terradish_direct_control_defaults <- function()
{
  list(
    factorization = "auto",
    solve_backend = "matrix",
    supernodal_min_vertices = 50000L,
    supernodal_max_rhs = 64L,
    perm = TRUE
  )
}

.terradish_resolve_direct_factorization <- function(control, n_vertices, n_rhs)
{
  factorization <- if (is.null(control$factorization)) "auto" else control$factorization
  factorization <- match.arg(factorization, c("auto", "simplicial_ldl", "simplicial_ll", "supernodal_ll"))

  if (!identical(factorization, "auto"))
    return(factorization)

  if (n_vertices >= as.integer(control$supernodal_min_vertices) &&
      n_rhs <= as.integer(control$supernodal_max_rhs))
    return("supernodal_ll")

  "simplicial_ldl"
}

.terradish_resolve_direct_solve_backend <- function(control)
{
  solve_backend <- if (is.null(control$solve_backend)) "matrix" else control$solve_backend
  match.arg(solve_backend, c("matrix", "cholmod_cpp", "cholmod_cpp_cached"))
}

.terradish_direct_signature <- function(control, factorization, solve_backend)
{
  list(
    factorization = factorization,
    perm = isTRUE(control$perm),
    solve_backend = solve_backend
  )
}

.terradish_solver_warning_is_fatal <- function(message)
{
  grepl("matrix not positive definite", message, ignore.case = TRUE) ||
    grepl("CHOLMOD", message, ignore.case = TRUE)
}

.terradish_with_solver_warning_guard <- function(expr, context)
{
  withCallingHandlers(
    expr,
    warning = function(w)
    {
      message <- conditionMessage(w)
      if (.terradish_solver_warning_is_fatal(message))
      {
        stop(
          context,
          " ",
          message,
          " Current conductance values produced a non-positive-definite reduced Laplacian.",
          call. = FALSE
        )
      }
    }
  )
}

.choleski_template <- function(s, factorization)
{
  templates <- s$choleski_templates
  if (!is.null(templates[[factorization]]))
    return(templates[[factorization]])
  if (identical(factorization, "simplicial_ldl") && !is.null(s$choleski))
    return(s$choleski)
  NULL
}

.normalize_solver_control <- function(solver, solver_control)
{
  solver_control <- if (is.null(solver_control)) NULL else as.list(solver_control)

  defaults <- switch(solver,
                     direct = .terradish_direct_control_defaults(),
                     amg = list(tol = 1e-8,
                                maxit = 400L,
                                coarse_enough = 1000L,
                                npre = 1L,
                                npost = 1L,
                                sa_relax = 1,
                                aggr_eps_strong = 0.08,
                                estimate_spectral_radius = TRUE,
                                power_iters = 4L,
                                reuse_preconditioner = TRUE,
                                reuse_preconditioner_max_age = Inf),
                     pcg_jacobi = list(tol = 1e-8, maxit = 1000L),
                     pcg = list(tol = 1e-8, maxit = 1000L),
                     block_cg = list(tol = 1e-8, maxit = 1000L),
                     stop("Unknown solver: ", solver))

  if (identical(solver, "amg") && !is.null(solver_control))
  {
    if (is.null(solver_control$tol) && !is.null(solver_control$tol_final))
      solver_control$tol <- solver_control$tol_final
    if (is.null(solver_control$maxit) && !is.null(solver_control$maxit_final))
      solver_control$maxit <- solver_control$maxit_final
  }

  control <- if (is.null(solver_control))
    defaults
  else
    modifyList(defaults, solver_control)

  if (identical(solver, "amg"))
  {
    max_age <- control$reuse_preconditioner_max_age
    if (length(max_age) != 1L || is.na(max_age) || max_age < 0)
      stop("`solver_control$reuse_preconditioner_max_age` must be a nonnegative number or `Inf`")
  }

  control
}

.terradish_solver_setup <- function(s, conductance, solver, solver_control = NULL, solver_reuse_state = NULL)
{
  requested_solver <- match.arg(solver, c("direct", "auto", "amg", "pcg", "pcg_jacobi", "block_cg"))
  resolution <- .terradish_resolve_solver(s, requested_solver, length(conductance), solver_control = solver_control)
  solver <- resolution$type
  control <- .normalize_solver_control(solver, resolution$solver_control)

  if (solver == "direct")
  {
    Qn <- .graph_reduced_laplacian(s, conductance)
    factorization <- .terradish_resolve_direct_factorization(control, resolution$n_vertices, resolution$n_rhs)
    solve_backend <- .terradish_resolve_direct_solve_backend(control)
    signature <- .terradish_direct_signature(control, factorization, solve_backend)

    if (identical(solve_backend, "cholmod_cpp_cached"))
    {
      can_reuse <- !is.null(solver_reuse_state) &&
        identical(solver_reuse_state$type, "direct") &&
        identical(solver_reuse_state$signature, signature) &&
        !is.null(solver_reuse_state$handle)

      setup_start <- proc.time()[["elapsed"]]
      handle <- .terradish_with_solver_warning_guard(
        if (can_reuse)
        {
          cholmod_direct_update(solver_reuse_state$handle, Qn)
          solver_reuse_state$handle
        }
        else
        {
          cholmod_direct_create(Qn, factorization, isTRUE(control$perm))
        },
        context = "Failed to prepare the direct solver."
      )
      setup_time <- proc.time()[["elapsed"]] - setup_start

      return(list(type = "direct",
                  requested_type = resolution$requested_type,
                  auto_reason = resolution$reason,
                  n_vertices = resolution$n_vertices,
                  n_rhs = resolution$n_rhs,
                  factor = NULL,
                  handle = handle,
                  factorization = factorization,
                  solve_backend = solve_backend,
                  signature = signature,
                  reused_factor_template = can_reuse,
                  setup_time = setup_time,
                  control = control))
    }

    can_reuse <- !is.null(solver_reuse_state) &&
      identical(solver_reuse_state$type, "direct") &&
      identical(solver_reuse_state$signature, signature) &&
      !is.null(solver_reuse_state$factor)
    template <- if (can_reuse) solver_reuse_state$factor else .choleski_template(s, factorization)

    setup_start <- proc.time()[["elapsed"]]
    factor <- .terradish_with_solver_warning_guard(
      if (is.null(template))
      {
        switch(factorization,
               simplicial_ldl = Cholesky(Qn, LDL = TRUE, super = FALSE, perm = isTRUE(control$perm)),
               simplicial_ll = Cholesky(Qn, LDL = FALSE, super = FALSE, perm = isTRUE(control$perm)),
               supernodal_ll = Cholesky(Qn, LDL = FALSE, super = TRUE, perm = isTRUE(control$perm)),
               stop("Unknown direct factorization mode: ", factorization))
      }
      else
        update(template, Qn),
      context = "Failed to prepare the direct solver."
    )
    setup_time <- proc.time()[["elapsed"]] - setup_start

    return(list(type = "direct",
                requested_type = resolution$requested_type,
                auto_reason = resolution$reason,
                n_vertices = resolution$n_vertices,
                n_rhs = resolution$n_rhs,
                factor = factor,
                factorization = factorization,
                solve_backend = solve_backend,
                signature = signature,
                reused_factor_template = can_reuse,
                setup_time = setup_time,
                control = control))
  }

  if (solver == "amg")
  {
    edge_pairs <- .graph_edge_pairs(s, length(conductance))
    signature <- .terradish_amg_reuse_signature(control)
    reuse_age <- as.integer(if (is.null(solver_reuse_state$reuse_age)) 0L else solver_reuse_state$reuse_age)
    reuse_max_age <- control$reuse_preconditioner_max_age
    can_reuse <- isTRUE(control$reuse_preconditioner) &&
      !is.null(solver_reuse_state) &&
      identical(solver_reuse_state$type, "amg") &&
      identical(solver_reuse_state$signature, signature) &&
      identical(as.integer(solver_reuse_state$n_vertices), resolution$n_vertices) &&
      (is.infinite(reuse_max_age) || reuse_age < as.integer(reuse_max_age))

    if (can_reuse)
    {
      amg_reduced_laplacian_rebuild(
        solver_reuse_state$handle,
        conductance,
        edge_pairs
      )
      return(list(type = "amg",
                  requested_type = resolution$requested_type,
                  auto_reason = resolution$reason,
                  n_vertices = resolution$n_vertices,
                  n_rhs = resolution$n_rhs,
                  handle = solver_reuse_state$handle,
                  conductance = conductance,
                  edge_pairs = edge_pairs,
                  control = control,
                  signature = signature,
                  reused_preconditioner = TRUE,
                  reuse_age = reuse_age + 1L))
    }

    return(list(type = "amg",
                requested_type = resolution$requested_type,
                auto_reason = resolution$reason,
                n_vertices = resolution$n_vertices,
                n_rhs = resolution$n_rhs,
                handle = amg_reduced_laplacian_create(
                  conductance,
                  edge_pairs,
                  tol = control$tol,
                  maxit = as.integer(control$maxit),
                  coarse_enough = as.integer(control$coarse_enough),
                  npre = as.integer(control$npre),
                  npost = as.integer(control$npost),
                  sa_relax = control$sa_relax,
                  aggr_eps_strong = control$aggr_eps_strong,
                  estimate_spectral_radius = isTRUE(control$estimate_spectral_radius),
                  power_iters = as.integer(control$power_iters)
                ),
                conductance = conductance,
                edge_pairs = edge_pairs,
                control = control,
                signature = signature,
                reused_preconditioner = FALSE,
                reuse_age = 0L))
  }

  list(type = solver,
       requested_type = resolution$requested_type,
       auto_reason = resolution$reason,
       n_vertices = resolution$n_vertices,
       n_rhs = resolution$n_rhs,
       conductance = conductance,
       edge_pairs = .graph_edge_pairs(s, length(conductance)),
       control = control)
}

.terradish_solver_solve <- function(solver_state, rhs, warm_start = NULL)
{
  common_info <- list(requested_type = solver_state$requested_type,
                      auto_reason = solver_state$auto_reason,
                      n_vertices = solver_state$n_vertices,
                      n_rhs = solver_state$n_rhs)

  if (solver_state$type == "direct")
  {
    solve_backend <- if (is.null(solver_state$solve_backend))
      "matrix"
    else
      solver_state$solve_backend

    if (identical(solve_backend, "cholmod_cpp_cached"))
    {
      out <- cholmod_direct_solve(solver_state$handle, as.matrix(rhs))
      solution <- out$solution
      solve_time <- out$solve_time
    }
    else if (identical(solve_backend, "cholmod_cpp"))
    {
      out <- cholmod_factor_solve(solver_state$factor, as.matrix(rhs))
      solution <- out$solution
      solve_time <- out$solve_time
    }
    else
    {
      solve_start <- proc.time()[["elapsed"]]
      solution <- .terradish_with_solver_warning_guard(
        solve(solver_state$factor, rhs),
        context = "Failed to solve the reduced Laplacian."
      )
      solve_time <- proc.time()[["elapsed"]] - solve_start
    }

    return(list(solution = solution,
                info = c(list(type = "direct",
                              factorization = solver_state$factorization,
                              solve_backend = solve_backend,
                              reused_factor_template = isTRUE(solver_state$reused_factor_template),
                              setup_time = solver_state$setup_time,
                              solve_time = solve_time),
                         common_info),
                warm_start = as.matrix(solution)))
  }

  if (solver_state$type == "amg")
  {
    out <- amg_reduced_laplacian_solve(
      solver_state$handle,
      as.matrix(rhs),
      x0 = if (is.null(warm_start)) NULL else as.matrix(warm_start),
      tol = solver_state$control$tol,
      maxit = as.integer(solver_state$control$maxit)
    )
    if (!all(out$converged))
      stop("AMG solver failed to converge for ", sum(!out$converged), " RHS column(s)")
    return(list(solution = out$solution,
                info = c(list(type = "amg",
                              converged = out$converged,
                              iterations = out$iterations,
                              residual_norm = out$residual_norm,
                              target_tol = solver_state$control$tol,
                              target_maxit = solver_state$control$maxit,
                              adaptive_phase = solver_state$control$adaptive_phase,
                              reused_preconditioner = isTRUE(solver_state$reused_preconditioner),
                              reuse_age = as.integer(solver_state$reuse_age),
                              setup_time = out$setup_time,
                              solve_time = out$solve_time,
                              n_reduced = out$n_reduced,
                              nnz = out$nnz),
                         common_info),
                warm_start = out$solution))
  }

  rhs <- as.matrix(rhs)
  if (!is.null(warm_start))
    warm_start <- as.matrix(warm_start)
  solver_fun <- switch(solver_state$type,
                       pcg = pcg_reduced_laplacian_ic,
                       pcg_jacobi = pcg_reduced_laplacian,
                       block_cg = block_cg_reduced_laplacian,
                       stop("Unknown solver type: ", solver_state$type))
  out <- solver_fun(rhs,
                    solver_state$conductance,
                    solver_state$edge_pairs,
                    x0 = warm_start,
                    tol = solver_state$control$tol,
                    maxit = as.integer(solver_state$control$maxit))
  if (!all(out$converged))
    stop(toupper(solver_state$type), " solver failed to converge for ", sum(!out$converged), " RHS column(s)")
  list(solution = out$solution,
       info = c(list(type = solver_state$type,
                     converged = out$converged,
                     iterations = out$iterations,
                     residual_norm = out$residual_norm,
                     target_tol = solver_state$control$tol,
                     target_maxit = solver_state$control$maxit),
                common_info),
       warm_start = out$solution)
}

.terradish_algorithm_derivative_results <- function(idx, state, cores, worker_libpaths = .libPaths())
{
  n_workers <- min(as.integer(cores), length(idx))
  splits <- split(idx, cut(idx, breaks = n_workers, labels = FALSE))
  cl <- makeCluster(length(splits))
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(cl, varlist = c("worker_libpaths"), envir = environment())
  clusterEvalQ(cl, {
    .libPaths(worker_libpaths)
    library(terradish)
    NULL
  })

  chunks <- parLapply(cl, splits, .terradish_algorithm_derivative_chunk, state = state)
  out <- unlist(chunks, recursive = FALSE)
  attr(out, "solver_info") <- lapply(chunks, attr, "solver_info")
  out
}

.terradish_solver_info_list <- function(info)
{
  if (is.null(info))
    return(list())
  if (!is.null(names(info)) && "type" %in% names(info))
    return(list(info))
  info
}

.terradish_solver_info_value <- function(info, name)
{
  value <- info[[name]]
  if (is.null(value))
    return(0)
  sum(as.numeric(value), na.rm = TRUE)
}

#' Likelihood of parameterized conductance surface
#'
#' Calculates likelihood, gradient, hessian, and partial derivatives of a
#' parameterized conductance surface, given a function mapping spatial data to
#' conductance and a function mapping resistance distance (covariance) to
#' genetic distance; using the algorithm in Pope (in prep).
#'
#' @param f A function of class 'conductance_model'
#' @param g A function of class 'measurement_model'
#' @param s An object of class \code{"terradish_graph"}
#' @param S A matrix of observed genetic distances
#' @param theta Parameters for conductance surface (e.g. inputs to 'f')
#' @param phi Optional warm-start values for the nuisance-parameter subproblem.
#' @param nu Effective Wishart degrees of freedom, passed through to
#'   measurement models that require it.
#' @param objective Compute negative loglikelihood?
#' @param gradient Compute gradient of negative loglikelihood wrt theta?
#' @param hessian Compute Hessian matrix of negative loglikelihood wrt theta?
#' @param partial Compute partial derivatives of negative loglikelihood wrt theta and spatial covariates/observed genetic distances
#' @param nonnegative Force regression-like 'measurement_model' to have nonnegative slope?
#' @param validate Numerical validation via 'numDeriv' (very slow, use for debugging small examples)
#' @param cores Number of worker processes to use for per-parameter derivative calculations. \code{1} evaluates serially.
#' @param curvature Which curvature to return in \code{hessian}. \code{"exact"}
#'   (default) returns the exact Hessian of the negative log-likelihood.
#'   \code{"gauss_newton"} returns the Gauss-Newton / Fisher-information
#'   approximation, obtained by dropping the two residual-weighted
#'   second-derivative terms (the second derivative of the resistance covariance
#'   \eqn{E} and of conductance with respect to \code{theta}). The Gauss-Newton
#'   curvature equals the Fisher information at the optimum, where it is positive
#'   semidefinite and yields a well-defined \code{vcov} (away from the optimum
#'   the measurement model's observed curvature can be indefinite). It equals the
#'   exact Hessian at a well-fitting optimum, and the gap between them measures
#'   model misspecification. It also requires only first derivatives of the
#'   conductance model, which is useful for Gaussian scale-aware and spline
#'   conductance models whose second derivatives are expensive or unstable. The
#'   number of linear solves is the same as for the exact Hessian.
#' @param solver Linear-system solver used for the reduced Laplacian. \code{"direct"} uses the cached sparse Cholesky factorization, \code{"auto"} conservatively chooses between the direct and AMG backends based on graph size and right-hand-side count, \code{"amg"} uses smoothed-aggregation algebraic multigrid preconditioned conjugate gradients, \code{"pcg"} uses incomplete-Cholesky preconditioned conjugate gradients, and \code{"pcg_jacobi"} keeps the older Jacobi-preconditioned prototype.
#' @param solver_control Optional named list of solver settings. For \code{solver = "direct"}, supported entries include \code{factorization} (\code{"auto"}, \code{"simplicial_ldl"}, \code{"simplicial_ll"}, or \code{"supernodal_ll"}), \code{solve_backend} (\code{"matrix"} or the experimental \code{"cholmod_cpp"} and \code{"cholmod_cpp_cached"} backends), \code{supernodal_min_vertices}, \code{supernodal_max_rhs}, and \code{perm}. For \code{solver = "auto"}, supported selection entries include \code{auto_direct_max_vertices}, \code{auto_amg_min_vertices}, and \code{auto_direct_max_rhs}. For \code{solver = "amg"}, supported entries include \code{tol}, \code{maxit}, \code{coarse_enough}, \code{npre}, \code{npost}, \code{sa_relax}, \code{aggr_eps_strong}, \code{estimate_spectral_radius}, \code{power_iters}, and \code{reuse_preconditioner}. For \code{solver = "pcg"} or \code{"pcg_jacobi"}, supported entries are \code{tol} and \code{maxit}.
#'   \code{reuse_preconditioner_max_age} can be set to a finite nonnegative
#'   value to periodically rebuild the AMG hierarchy instead of reusing it
#'   indefinitely.
#'   Direct supernodal factorizations can benefit from a threaded BLAS, but the
#'   relevant thread counts are controlled by the external R/BLAS build rather
#'   than by \code{terradish_algorithm()}.
#'   If the direct solver stops with a CHOLMOD or
#'   non-positive-definite reduced Laplacian message, the current conductance
#'   values made the graph numerically singular. Check for missing or infinite
#'   raster covariates, use scaled covariates and moderate starting values near
#'   zero, and try \code{solver = "auto"} or \code{solver = "amg"} while
#'   diagnosing the fit.
#' @param solver_warm_start Optional initial guess for the reduced-system solve. This is primarily useful for iterative solvers when evaluating nearby parameter values.
#' @param solver_reuse_state Optional reusable solver state returned by a prior
#'   \code{terradish_algorithm()} call. This is used to reuse AMG hierarchy
#'   information or compatible direct CHOLMOD factorization state across nearby
#'   evaluations.
#'
#' @return A list containing at a minimum:
#'  \item{covariance}{rows/columns of the generalized inverse of the graph Laplacian for a subset of target vertices}
#' Additionally, if 'objective == TRUE':
#'  \item{objective}{(if 'objective') the negative loglikelihood}
#'  \item{phi}{(if 'objective') fitted values of the nuisance parameters of 'g'}
#'  \item{boundary}{(if 'objective') is the solution on the boundary (e.g. no genetic structure)?}
#'  \item{fitted}{(if 'objective') matrix of expected genetic distances among target vertices}
#'  \item{gradient}{(if 'gradient') gradient of negative loglikelihood with respect to theta}
#'  \item{hessian}{(if 'hessian') Hessian matrix of the negative loglikelihood with respect to theta}
#'  \item{partial_X}{(if 'partial') Jacobian of the gradient with respect to the spatial covariates}
#'  \item{partial_S}{(if 'partial') Jacobian of the gradient with respect to the observed genetic distances}
#'
#' @examples
#' library(terra)
#'
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.forestcover <- terra::unwrap(melip.forestcover)
#' melip.coords <- terra::unwrap(melip.coords)
#'
#' covariates <- c(melip.altitude, melip.forestcover)
#' names(covariates) <- c("altitude", "forestcover")
#' surface <- conductance_surface(covariates, melip.coords, directions = 8)
#' conductance_model <- loglinear_conductance(~ altitude + forestcover, surface$x)
#'
#' terradish_algorithm(conductance_model, terradish::leastsquares,
#'                  surface, ifelse(melip.Fst < 0, 0, melip.Fst),
#'                  nu = 1000, theta = c(-0.3, 0.3))
#'
#' @export
terradish_algorithm <- function(f, g, s, S, theta, nu = NULL, phi = NULL, objective = TRUE, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE, validate = FALSE, cores = 1L, curvature = c("exact", "gauss_newton"), solver = c("direct", "auto", "amg", "pcg", "pcg_jacobi", "block_cg"), solver_control = NULL, solver_warm_start = NULL, solver_reuse_state = NULL)
{
  stopifnot(inherits(f, c("terradish_conductance_model",
                          "radish_conductance_model")))
  stopifnot(inherits(g, c("terradish_measurement_model",
                          "radish_measurement_model")))
  stopifnot(inherits(s, c("terradish_graph", "radish_graph")))

  stopifnot(is.matrix(S)     )
  stopifnot(is.numeric(theta))

  stopifnot(length(s$demes) == nrow(S)  )
  stopifnot(        nrow(S) == ncol(S)  )
  solver <- match.arg(solver)
  curvature <- match.arg(curvature)
  gauss_newton <- identical(curvature, "gauss_newton")

  symm <- function(X) (X + t(X))/2

  # conductance
  C <- f(theta)
  df__dtheta_matrix <- .conductance_df_dtheta_matrix(C, theta)

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- length(C$conductance)
  reduced_index <- .graph_reduced_index(s, N)
  Zn   <- .graph_rhs(s, N)
  solver_state <- .terradish_solver_setup(s, C$conductance, solver = solver, solver_control = solver_control, solver_reuse_state = solver_reuse_state)
  solve_main <- .terradish_solver_solve(solver_state, Zn, warm_start = solver_warm_start)
  derivative_solver_info <- list()
  G    <- solve_main$solution
  Gmat <- as.matrix(G)
  tG   <- t(Gmat)
  E    <- graph_rhs_crossprod(s$demes, N, Gmat)

  if (objective || gradient || hessian)
  {
    # measurement model
    E_dense <- as.matrix(E)
    subproblem <- radish_subproblem(g = g, E = E_dense, S = S, nu = nu, phi = phi,
                                    nonnegative = nonnegative,
                                    control = NewtonRaphsonControl(verbose = FALSE, 
                                                                   ftol = 1e-10, 
                                                                   ctol = 1e-10))
    phi        <- subproblem$phi
    loglik     <- subproblem$loglikelihood

    # gradient calculation
    grad      <- rep(0, length(theta))
    hess      <- matrix(0, length(theta), length(theta))
    partial_X <- NULL
    partial_S <- NULL
    if (gradient || hessian || partial)
    {
      dl_dE    <- subproblem$gradient 
      dl_dQnG  <- dl_dE %*% tG
      dl_dC    <- backpropagate_laplacian_to_conductance(dl_dQnG, tG, s$adj)
      grad <- c(crossprod(df__dtheta_matrix, c(dl_dC)))

      # hessian and mixed partial derivative calculations
      if (hessian || partial)
      {
        if (partial)
        {
          partial_X <- array(0, c(N, length(theta), length(theta)))
          partial_S <- array(0, c(nrow(S), ncol(S), length(theta)))
        }
        idx <- seq_along(theta)
        can_parallel <- as.integer(cores) > 1L &&
          length(idx) > 1L &&
          !identical(solver_state$type, "amg") &&
          identical(environmentName(environment(terradish_algorithm)), "namespace:terradish")

        if (can_parallel)
        {
          derivative_state <- list(C = C,
                                   s = s,
                                   G = Gmat,
                                   tG = tG,
                                   Zn = Zn,
                                   solver_state = solver_state,
                                   df__dtheta_matrix = df__dtheta_matrix,
                                   subproblem = subproblem,
                                   dl_dE = dl_dE,
                                   dl_dC = dl_dC,
                                   theta = theta,
                                   partial = partial,
                                   gauss_newton = gauss_newton,
                                   N = N)
          derivative_results <- .terradish_algorithm_derivative_results(
            idx = idx,
            state = derivative_state,
            cores = cores,
            worker_libpaths = .libPaths()
          )
        }
        else
        {
          derivative_state <- list(C = C,
                                   s = s,
                                   G = Gmat,
                                   tG = tG,
                                   Zn = Zn,
                                   solver_state = solver_state,
                                   df__dtheta_matrix = df__dtheta_matrix,
                                   subproblem = subproblem,
                                   dl_dE = dl_dE,
                                   dl_dC = dl_dC,
                                   theta = theta,
                                   partial = partial,
                                   gauss_newton = gauss_newton,
                                   N = N)
    derivative_results <- .terradish_algorithm_derivative_chunk(idx, derivative_state)
        }
        derivative_solver_info <- .terradish_solver_info_list(
          attr(derivative_results, "solver_info")
        )

        for (res in derivative_results)
        {
          k <- res$k
          hess[k, ] <- res$hess_row
          if (partial)
          {
            partial_X[, k, ] <- res$partial_X_k
            partial_S[, , k] <- res$partial_S_k
          }
        }
        if (gauss_newton)
          hess <- symm(hess)
      }
    }
  }

  # numerical validation
  if (validate)
  {
    num_gradient <- .numderiv_grad(function(theta) 
                                   terradish_algorithm(f = f,
                                                    g = g, 
                                                    s = s,
                                                    S = S, 
                                                    theta = theta)$objective, 
                                   theta, method = "Richardson")

    num_hessian  <- .numderiv_hessian(function(theta) 
                                      terradish_algorithm(f = f,
                                                       g = g, 
                                                      s = s,
                                                      S = S, 
                                                      theta = theta)$objective, 
                                      theta, method = "Richardson")

    warning("Numerical validation of `partial_X` is disabled for formula-based conductance models.",
            call. = FALSE)
    num_partial_X <- NULL

    num_partial_S <- array(t(.numderiv_jacobian(function(S) 
                                                terradish_algorithm(f = f,
                                                                 g = g, 
                                                                 s = s,
                                                                 S = S, 
                                                                 theta = theta)$gradient, 
                                                S, method = "simple")), 
                           c(nrow(S), ncol(S), length(theta)))
  }

  solver_infos <- c(list(solve_main$info), derivative_solver_info)
  algorithm_diagnostics <- list(
    solver_setups = 1L,
    solver_solves = length(solver_infos),
    solver_setup_time = .terradish_solver_info_value(solve_main$info, "setup_time"),
    solver_solve_time = sum(vapply(solver_infos,
                                    .terradish_solver_info_value,
                                    numeric(1),
                                    name = "solve_time"))
  )

  list (covariance    = E,
         objective     = if(!objective) NULL else loglik,
         phi           = if(!objective) NULL else phi,
         phi_hessian   = if(!objective) NULL else subproblem$fit$hessian,
         boundary      = if(!objective) NULL else subproblem$boundary, # the solution is on the boundary (e.g. no genetic structure) so all derivatives wrt theta are 0
         fitted        = if(!objective) NULL else subproblem$fit$fitted,
         gradient      = if(!gradient)  NULL else grad * (1 - subproblem$boundary), # wrt theta
        hessian       = if(!hessian)   NULL else hess * (1 - subproblem$boundary), # wrt theta
        partial_X     = if(!partial)   NULL else partial_X * (1 - subproblem$boundary), # partial_X[i,l,k] is \frac{\partial^2 L(theta,x)}{\partial theta_l \partial x_{ik}}
         partial_S     = if(!partial)   NULL else partial_S * (1 - subproblem$boundary), # partial_S[i,j,k] is \frac{\partial^2 L(theta,x)}{\partial theta_k \partial S_{ij}}
         num_gradient  = if(!validate)  NULL else num_gradient,
         num_hessian   = if(!validate)  NULL else num_hessian,
         num_partial_X = if(!validate)  NULL else num_partial_X,
         num_partial_S = if(!validate)  NULL else num_partial_S,
         solver_info   = solve_main$info,
         algorithm_diagnostics = algorithm_diagnostics,
         solver_warm_start = solve_main$warm_start,
      solver_reuse_state = switch(
            solver_state$type,
            amg = list(
              type = "amg",
             handle = solver_state$handle,
             signature = solver_state$signature,
             reuse_age = solver_state$reuse_age,
             n_vertices = solver_state$n_vertices
           ),
           direct = if (identical(solver_state$solve_backend, "cholmod_cpp_cached"))
           {
             list(
               type = "direct",
               handle = solver_state$handle,
               signature = solver_state$signature
             )
           }
           else
           {
             list(
               type = "direct",
               factor = solver_state$factor,
               signature = solver_state$signature
             )
           },
            NULL
          ))
}

#' Legacy radish algorithm wrapper
#'
#' Deprecated compatibility wrapper retained for legacy code that still calls
#' \code{radish_algorithm()} directly.
#'
#' @param ... Arguments passed through the deprecated
#'   \code{radish_algorithm()} compatibility wrapper to
#'   \code{\link{terradish_algorithm}}.
#' @name legacy_radish_algorithm_wrapper
#' @keywords internal
NULL

#' @rdname legacy_radish_algorithm_wrapper
#' @export
radish_algorithm <- function(...)
{
  .terradish_deprecate("radish_algorithm", "terradish_algorithm")
  .terradish_forward_call(match.call(), "terradish_algorithm")
}

# Internal compatibility aliases retained for housekeeping-sized refactors.
.radish_algorithm_derivative_chunk <- .terradish_algorithm_derivative_chunk
.radish_auto_solver_defaults <- .terradish_auto_solver_defaults
.radish_resolve_solver <- .terradish_resolve_solver
.radish_amg_reuse_signature <- .terradish_amg_reuse_signature
.radish_direct_control_defaults <- .terradish_direct_control_defaults
.radish_resolve_direct_factorization <- .terradish_resolve_direct_factorization
.radish_direct_signature <- .terradish_direct_signature
.radish_solver_setup <- .terradish_solver_setup
.radish_solver_solve <- .terradish_solver_solve
.radish_algorithm_derivative_results <- .terradish_algorithm_derivative_results
