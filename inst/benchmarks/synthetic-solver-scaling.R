#!/usr/bin/env Rscript

# Synthetic solver-scaling benchmark for terradish.
#
# This script uses:
# - multiScaleR::sim_rast() to generate terra-native covariate rasters
# - PopGenReport to simulate population-genetic data and derive pairwise FST
# - terradish to fit resistance models with different solver backends
#
# The main goal is to identify when the exact sparse-Cholesky solver stops being
# competitive relative to iterative methods as graph size increases.

CONFIG <- list(
  raster_dims = c(100L, 250L, 500L),
  n_pops_values = c(25L, 50L),
  n_candidate_pops = 50L,
  resolution = 1,
  autocorr_range1 = 8,
  autocorr_range2 = 25,
  sill = 10,
  n_ind = 20L,
  sex_ratio = 0.5,
  n_loci = 25L,
  n_allels = 12L,
  steps = 150L,
  n_offspring = 2L,
  mig_rate = 0.10,
  disp_quantile = 0.10,
  mutation_rate = 0,
  directions = 8L,
  formula = fst ~ cont1 + cont2,
  measurement_model = "generalized_wishart",
  nu = NULL,
  true_theta = c(cont1 = 0.8, cont2 = -0.8),
  reps = 2L,
  max_sim_attempts = 5L,
  timeout_sec = 600,
  terradish_cores = 1L,
  solvers = c("direct", "auto", "amg"),
  auto_control = list(
    auto_direct_max_vertices = 750000L,
    auto_amg_min_vertices = 1500000L,
    auto_direct_max_rhs = 64L,
    adaptive = TRUE,
    tol_early = 1e-4,
    tol_mid = 1e-6,
    tol_final = 1e-8,
    maxit_early = 100L,
    maxit_mid = 250L,
    maxit_final = 400L,
    warmup_evals = 3L,
    coarse_enough = 1000L,
    npre = 1L,
    npost = 1L,
    sa_relax = 1,
    aggr_eps_strong = 0.08,
    estimate_spectral_radius = TRUE,
    power_iters = 4L,
    reuse_preconditioner = TRUE
  ),
  amg_control = list(
    adaptive = TRUE,
    tol_early = 1e-4,
    tol_mid = 1e-6,
    tol_final = 1e-8,
    maxit_early = 100L,
    maxit_mid = 250L,
    maxit_final = 400L,
    warmup_evals = 3L,
    coarse_enough = 1000L,
    npre = 1L,
    npost = 1L,
    sa_relax = 1,
    aggr_eps_strong = 0.08,
    estimate_spectral_radius = TRUE,
    power_iters = 4L,
    reuse_preconditioner = TRUE
  ),
  pcg_control = list(tol = 1e-8, maxit = 5000L),
  seed = 1L
)

with_elapsed_timeout <- function(expr, timeout_sec)
{
  old_limit <- base::getOption("expressions")
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  force(expr)
}

sample_focal_coords <- function(r, n_pops, seed)
{
  set.seed(seed)
  cells <- terra::spatSample(r[[1]], size = n_pops, method = "random",
                             na.rm = TRUE, cells = TRUE, xy = TRUE,
                             values = FALSE, as.points = FALSE)
  coords <- as.data.frame(cells[, c("x", "y"), drop = FALSE])
  rownames(coords) <- sprintf("pop_%02d", seq_len(nrow(coords)))
  coords
}

candidate_population_count <- function(config)
{
  max(as.integer(config$n_candidate_pops), 2L * as.integer(config$n_pops))
}

build_cost_matrix <- function(coords)
{
  as.matrix(stats::dist(coords[, c("x", "y"), drop = FALSE]))
}

prepare_covariates <- function(rast)
{
  covariates <- c(terra::scale(rast[["cont1"]]), terra::scale(rast[["cont2"]]))
  names(covariates) <- c("cont1", "cont2")
  covariates
}

true_cost_matrix <- function(surface, config)
{
  theta <- matrix(unname(config$true_theta), nrow = 1)
  colnames(theta) <- names(config$true_theta)
  cost <- terradish::terradish_distance(
    theta = theta,
    formula = ~ cont1 + cont2,
    data = surface,
    conductance_model = terradish::loglinear_conductance,
    covariance = FALSE,
    cores = 1L
  )
  as.matrix(cost$distance[, , 1])
}

simulate_fst <- function(coords, cost.mat, config)
{
  simpops <- PopGenReport::init.popgensim(
    n.pops = nrow(coords),
    n.ind = config$n_ind,
    sex.ratio = config$sex_ratio,
    n.loci = config$n_loci,
    n.allels = config$n_allels,
    locs = coords,
    n.cov = 3
  )

  disp.max <- as.numeric(stats::quantile(cost.mat[upper.tri(cost.mat)],
                                         probs = config$disp_quantile,
                                         na.rm = TRUE))

  pops <- PopGenReport::run.popgensim(
    simpops = simpops,
    steps = config$steps,
    cost.mat = cost.mat,
    n.offspring = config$n_offspring,
    n.ind = config$n_ind,
    mig.rate = config$mig_rate,
    disp.max = disp.max,
    disp.rate = config$disp_quantile,
    n.allels = config$n_allels,
    mut.rate = config$mutation_rate,
    n.cov = 3
  )

  gi <- PopGenReport::pops2genind(pops)
  fst <- PopGenReport::pairwise.fstb(gi)
  fst[is.na(fst)] <- 0
  fst
}

simulate_surface <- function(dim, config)
{
  rast <- multiScaleR::sim_rast(
    dim = dim,
    resolution = config$resolution,
    autocorr_range1 = config$autocorr_range1,
    autocorr_range2 = config$autocorr_range2,
    sill = config$sill,
    plot = FALSE,
    user_seed = config$seed + dim
  )

  names(rast) <- c("bin1", "bin2", "cont1", "cont2")
  covariates <- prepare_covariates(rast)
  for (attempt in seq_len(config$max_sim_attempts))
  {
    coords <- sample_focal_coords(rast, candidate_population_count(config),
                                  seed = config$seed + dim + 1000L * attempt)
    surface_candidate <- terradish::conductance_surface(covariates,
                                                        coords,
                                                        directions = config$directions)
    cost.mat <- true_cost_matrix(surface_candidate, config)
    fst <- simulate_fst(coords, cost.mat, config)
    keep <- rownames(fst)
    if (!is.null(keep) && length(keep))
    {
      if (all(grepl("^[0-9]+$", keep)))
      {
        coords <- coords[as.integer(keep), , drop = FALSE]
        cost.mat <- cost.mat[as.integer(keep), as.integer(keep), drop = FALSE]
      }
      else
      {
        coords <- coords[keep, , drop = FALSE]
        cost.mat <- cost.mat[keep, keep, drop = FALSE]
      }
    }
    if (nrow(coords) >= config$n_pops)
      break
  }
  if (nrow(coords) != nrow(fst))
    stop("Simulated focal coordinates and FST matrix are misaligned")
  if (nrow(coords) < config$n_pops)
    stop("Genetic simulation returned fewer than the requested number of focal populations")
  selected <- sort(sample(seq_len(nrow(coords)), size = config$n_pops, replace = FALSE))
  coords <- coords[selected, , drop = FALSE]
  fst <- fst[selected, selected, drop = FALSE]
  surface <- terradish::conductance_surface(covariates,
                                            coords,
                                            directions = config$directions)

  list(rast = rast,
       covariates = covariates,
       coords = coords,
       fst = fst,
       surface = surface)
}

fit_once <- function(surface_case, solver, config)
{
  fit <- NULL
  error_message <- NULL
  measurement_model <- get(config$measurement_model, envir = asNamespace("terradish"))
  nu <- if (is.null(config$nu)) config$n_loci else config$nu
  solver_control <- switch(
    solver,
    auto = config$auto_control,
    amg = config$amg_control,
    pcg = config$pcg_control,
    pcg_jacobi = config$pcg_control,
    NULL
  )
  elapsed <- tryCatch(
    {
      fst <- surface_case$fst
      fit_formula <- config$formula
      started <- proc.time()[["elapsed"]]
      fit <- with_elapsed_timeout(
        terradish::terradish(
          fit_formula,
          data = surface_case$surface,
          conductance_model = terradish::loglinear_conductance,
          measurement_model = measurement_model,
          nu = nu,
          control = terradish::NewtonRaphsonControl(maxit = 10, verbose = FALSE),
          cores = config$terradish_cores,
          solver = solver,
          solver_control = solver_control
        ),
        timeout_sec = config$timeout_sec
      )
      proc.time()[["elapsed"]] - started
    },
    error = function(e)
    {
      error_message <<- conditionMessage(e)
      NA_real_
    }
  )

  if (is.null(fit))
  {
    return(list(
      elapsed = elapsed,
      timed_out = is.na(elapsed),
      error = error_message,
      logLik = NA_real_,
      df = NA_real_,
      boundary = NA,
      solver_info = NULL
    ))
  }

  list(
    elapsed = elapsed,
    timed_out = FALSE,
    error = NA_character_,
    logLik = fit$loglik,
    df = fit$df,
    boundary = isTRUE(fit$fit$boundary),
    solver_info = fit$fit$solver_info
  )
}

benchmark_case <- function(dim, n_pops, config)
{
  cat(sprintf("\nPreparing synthetic case for dim = %d, n_pops = %d\n", dim, n_pops))
  config$n_pops <- as.integer(n_pops)
  surface_case <- tryCatch(simulate_surface(dim, config), error = identity)
  if (inherits(surface_case, "error"))
  {
    return(do.call(rbind, lapply(config$solvers, function(solver)
      data.frame(
        dim = dim,
        n_vertices = NA_integer_,
        n_pops = n_pops,
        solver = solver,
        rep = NA_integer_,
        elapsed_sec = NA_real_,
        timed_out = NA,
        logLik = NA_real_,
        df = NA_real_,
        boundary = NA,
        error = conditionMessage(surface_case),
        max_solver_iterations = NA_integer_,
        requested_solver = solver,
        resolved_solver = NA_character_,
        auto_reason = "",
        reused_preconditioner = NA,
        reuse_age = NA_integer_,
        solver_setup_time = NA_real_,
        solver_solve_time = NA_real_,
        stringsAsFactors = FALSE
      ))))
  }
  n_vertices <- nrow(surface_case$surface$x)

  results <- vector("list", length(config$solvers) * config$reps)
  idx <- 1L
  for (solver in config$solvers)
  {
    for (rep in seq_len(config$reps))
    {
      cat(sprintf("  solver = %s, rep = %d\n", solver, rep))
      run <- fit_once(surface_case, solver, config)
      results[[idx]] <- data.frame(
        dim = dim,
        n_vertices = n_vertices,
        n_pops = nrow(surface_case$fst),
        solver = solver,
        rep = rep,
        elapsed_sec = run$elapsed,
        timed_out = run$timed_out,
        logLik = run$logLik,
        df = run$df,
        boundary = isTRUE(run$boundary),
        error = if (is.null(run$error) || is.na(run$error)) "" else run$error,
        max_solver_iterations = if (!is.null(run$solver_info$iterations))
          max(run$solver_info$iterations) else NA_integer_,
        requested_solver = if (!is.null(run$solver_info$requested_type))
          run$solver_info$requested_type else solver,
        resolved_solver = if (!is.null(run$solver_info$type))
          run$solver_info$type else NA_character_,
        auto_reason = if (!is.null(run$solver_info$auto_reason))
          run$solver_info$auto_reason else "",
        reused_preconditioner = if (!is.null(run$solver_info$reused_preconditioner))
          isTRUE(run$solver_info$reused_preconditioner) else NA,
        reuse_age = if (!is.null(run$solver_info$reuse_age))
          as.integer(run$solver_info$reuse_age) else NA_integer_,
        solver_setup_time = if (!is.null(run$solver_info$setup_time))
          run$solver_info$setup_time else NA_real_,
        solver_solve_time = if (!is.null(run$solver_info$solve_time))
          run$solver_info$solve_time else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, results)
}

summarize_results <- function(results)
{
  split_results <- split(results, list(results$dim, results$n_pops, results$solver), drop = TRUE)
  out <- lapply(split_results, function(x)
  {
    data.frame(
      dim = x$dim[[1]],
      n_vertices = x$n_vertices[[1]],
      n_pops = x$n_pops[[1]],
      solver = x$solver[[1]],
      n_runs = nrow(x),
      n_timeouts = sum(x$timed_out, na.rm = TRUE),
      n_boundary = sum(x$boundary, na.rm = TRUE),
      median_elapsed_sec = stats::median(x$elapsed_sec, na.rm = TRUE),
      mean_elapsed_sec = mean(x$elapsed_sec, na.rm = TRUE),
      median_setup_sec = stats::median(x$solver_setup_time, na.rm = TRUE),
      median_solve_sec = stats::median(x$solver_solve_time, na.rm = TRUE),
      max_solver_iterations = if (all(is.na(x$max_solver_iterations))) NA_integer_
        else max(x$max_solver_iterations, na.rm = TRUE),
      resolved_solver = paste(unique(stats::na.omit(x$resolved_solver)), collapse = ","),
      any_reused_preconditioner = any(x$reused_preconditioner, na.rm = TRUE),
      max_reuse_age = if (all(is.na(x$reuse_age))) NA_integer_
        else max(x$reuse_age, na.rm = TRUE),
      logLik = if (all(is.na(x$logLik))) NA_real_ else x$logLik[which.max(!is.na(x$logLik))],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

main <- function(config = CONFIG)
{
  if (!requireNamespace("terradish", quietly = TRUE))
    stop("Package `terradish` must be installed.")
  if (!requireNamespace("multiScaleR", quietly = TRUE))
    stop("Package `multiScaleR` must be installed.")
  if (!requireNamespace("PopGenReport", quietly = TRUE))
    stop("Package `PopGenReport` must be installed.")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package `terra` must be installed.")

  if (any(config$n_pops_values < 25L))
    stop("Use at least 25 populations for this benchmark design.")
  if (config$n_candidate_pops < max(config$n_pops_values))
    stop("`n_candidate_pops` must be at least as large as the requested population counts.")

  cases <- expand.grid(dim = config$raster_dims,
                       n_pops = config$n_pops_values,
                       KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
  all_results <- lapply(seq_len(nrow(cases)), function(i)
    benchmark_case(dim = cases$dim[[i]], n_pops = cases$n_pops[[i]], config = config))
  all_results <- do.call(rbind, all_results)
  summary <- summarize_results(all_results)

  cat("\nRaw results:\n")
  print(all_results, row.names = FALSE)

  cat("\nSummary:\n")
  print(summary, row.names = FALSE)

  invisible(list(raw = all_results, summary = summary))
}

if (sys.nframe() == 0L)
  main()
