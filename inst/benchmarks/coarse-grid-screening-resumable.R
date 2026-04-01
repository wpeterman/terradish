#!/usr/bin/env Rscript

# Resumable benchmark for coarse-raster grid screening in terradish.
#
# This script checkpoints:
# - the simulated synthetic landscape
# - the exact grid result
# - each coarse-raster grid factor
# - optional exact-fit runs initialized from the best coarse-grid point
#
# The goal is to make larger runs (e.g. 1250x1250, 1500x1500) restartable
# without losing finished work if the session times out.

CONFIG <- list(
  dim = 1250L,
  n_pops = 25L,
  n_candidate_pops = 50L,
  factors = c(2L, 3L, 4L),
  center_mode = "true_theta",   # one of "true_theta", "exact_fit"
  theta_half_width = 0.25,
  grid_points = 3L,
  benchmark_starts = TRUE,
  output_dir = "C:/temp/terradish-coarse-screening-1250",
  seed = 23L
)

.synthetic_helpers <- new.env(parent = globalenv())
source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/inst/benchmarks/synthetic-solver-scaling.R",
       local = .synthetic_helpers)

.coarse_screening_defaults <- function()
{
  .synthetic_helpers$CONFIG
}

.merge_screening_config <- function(config)
{
  defaults <- .coarse_screening_defaults()
  defaults$seed <- config$seed
  defaults$n_pops <- as.integer(config$n_pops)
  defaults$n_candidate_pops <- as.integer(config$n_candidate_pops)
  defaults
}

.elapsed_run <- function(expr)
{
  gc()
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(elapsed = unname(proc.time()[["elapsed"]] - started),
       value = value)
}

.ensure_output_dir <- function(path)
{
  if (!dir.exists(path))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.paths <- function(output_dir)
{
  list(
    surface_case = file.path(output_dir, "surface_case.rds"),
    surface_raster = file.path(output_dir, "surface_raster_covariates.rds"),
    surface_selected = file.path(output_dir, "surface_selected.rds"),
    exact_fit = file.path(output_dir, "exact_fit.rds"),
    exact_grid = file.path(output_dir, "exact_grid.rds"),
    grid_raw = file.path(output_dir, "grid_results.csv"),
    grid_summary = file.path(output_dir, "grid_summary.csv"),
    start_raw = file.path(output_dir, "start_results.csv"),
    start_summary = file.path(output_dir, "start_summary.csv")
  )
}

.grid_rds_path <- function(output_dir, mode)
  file.path(output_dir, paste0("grid_", mode, ".rds"))

.surface_attempt_path <- function(output_dir, attempt)
  file.path(output_dir, sprintf("surface_attempt_%02d.rds", as.integer(attempt)))

.surface_coords_path <- function(output_dir, attempt)
  file.path(output_dir, sprintf("surface_attempt_%02d_coords.rds", as.integer(attempt)))

.surface_cost_path <- function(output_dir, attempt)
  file.path(output_dir, sprintf("surface_attempt_%02d_cost.rds", as.integer(attempt)))

.surface_fst_path <- function(output_dir, attempt)
  file.path(output_dir, sprintf("surface_attempt_%02d_fst.rds", as.integer(attempt)))

.read_or_empty <- function(path, prototype)
{
  if (file.exists(path))
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  prototype
}

.write_csv <- function(df, path)
{
  utils::write.csv(df, path, row.names = FALSE)
}

.pack_spat <- function(x)
{
  if (inherits(x, "PackedSpatRaster"))
    return(x)
  if (inherits(x, "SpatRaster"))
    return(terra::wrap(x))
  x
}

.unpack_spat <- function(x)
{
  if (inherits(x, "PackedSpatRaster"))
    return(terra::unwrap(x))
  x
}

.pack_surface_case <- function(surface_case)
{
  surface_case$rast <- .pack_spat(surface_case$rast)
  surface_case$covariates <- .pack_spat(surface_case$covariates)
  surface_case
}

.unpack_surface_case <- function(surface_case)
{
  surface_case$rast <- .unpack_spat(surface_case$rast)
  surface_case$covariates <- .unpack_spat(surface_case$covariates)
  surface_case
}

.build_theta_grid <- function(center, half_width, grid_points)
{
  theta <- as.matrix(expand.grid(
    cont2 = seq(center[["cont2"]] - half_width, center[["cont2"]] + half_width, length.out = grid_points),
    cont1 = seq(center[["cont1"]] - half_width, center[["cont1"]] + half_width, length.out = grid_points)
  ))
  theta[, c("cont1", "cont2"), drop = FALSE]
}

.grid_result_row <- function(label, run, surface_case, baseline = NULL)
{
  grid <- run$value
  data.frame(
    mode = label,
    elapsed_sec = run$elapsed,
    max_logLik = max(grid$loglik, na.rm = TRUE),
    min_logLik = min(grid$loglik, na.rm = TRUE),
    spearman = if (is.null(baseline)) 1 else suppressWarnings(stats::cor(baseline$loglik, grid$loglik, method = "spearman")),
    max_abs_logLik_diff = if (is.null(baseline)) 0 else max(abs(baseline$loglik - grid$loglik), na.rm = TRUE),
    mean_abs_logLik_diff = if (is.null(baseline)) 0 else mean(abs(baseline$loglik - grid$loglik), na.rm = TRUE),
    top_match = if (is.null(baseline)) TRUE else which.max(baseline$loglik) == which.max(grid$loglik),
    approx_type = grid$approximation$type,
    approx_used = isTRUE(grid$approximation$used),
    full_vertices = if (!is.null(grid$approximation$full_vertices)) grid$approximation$full_vertices else nrow(surface_case$surface$x),
    coarse_vertices = if (!is.null(grid$approximation$coarse_vertices)) grid$approximation$coarse_vertices else nrow(surface_case$surface$x),
    stringsAsFactors = FALSE
  )
}

.summarize_grid_results <- function(grid_results)
{
  out <- grid_results
  exact_elapsed <- out$elapsed_sec[out$mode == "exact"][1]
  out$speedup_vs_exact <- exact_elapsed / out$elapsed_sec
  out$vertex_fraction <- out$coarse_vertices / out$full_vertices
  out
}

.fit_result_row <- function(label, run, baseline)
{
  fit <- run$value
  data.frame(
    start = label,
    elapsed_sec = run$elapsed,
    logLik = unname(fit$loglik),
    max_abs_theta_diff = max(abs(fit$mle$theta - baseline$mle$theta)),
    stringsAsFactors = FALSE
  )
}

.summarize_fit_results <- function(fit_results)
{
  out <- fit_results
  baseline_elapsed <- out$elapsed_sec[out$start == "default"][1]
  out$speedup_vs_default <- baseline_elapsed / out$elapsed_sec
  out
}

.load_surface_raster <- function(paths, config, base_config)
{
  if (file.exists(paths$surface_raster))
  {
    raster_bundle <- readRDS(paths$surface_raster)
    raster_bundle$rast <- .unpack_spat(raster_bundle$rast)
    raster_bundle$covariates <- .unpack_spat(raster_bundle$covariates)
    return(raster_bundle)
  }

  cat(sprintf("Simulating %dx%d raster and covariates...\n", config$dim, config$dim))
  rast <- multiScaleR::sim_rast(
    dim = as.integer(config$dim),
    resolution = base_config$resolution,
    autocorr_range1 = base_config$autocorr_range1,
    autocorr_range2 = base_config$autocorr_range2,
    sill = base_config$sill,
    plot = FALSE,
    user_seed = base_config$seed + as.integer(config$dim)
  )
  names(rast) <- c("bin1", "bin2", "cont1", "cont2")
  raster_bundle <- list(
    rast = rast,
    covariates = .synthetic_helpers$prepare_covariates(rast)
  )
  saveRDS(list(
    rast = .pack_spat(raster_bundle$rast),
    covariates = .pack_spat(raster_bundle$covariates)
  ), paths$surface_raster)
  raster_bundle
}

.load_surface_attempt <- function(paths, output_dir, attempt, raster_bundle, config, base_config)
{
  attempt_path <- .surface_attempt_path(output_dir, attempt)
  if (file.exists(attempt_path))
    return(readRDS(attempt_path))

  coords_path <- .surface_coords_path(output_dir, attempt)
  cost_path <- .surface_cost_path(output_dir, attempt)
  fst_path <- .surface_fst_path(output_dir, attempt)

  if (file.exists(coords_path))
  {
    coords <- readRDS(coords_path)
  }
  else
  {
    cat(sprintf("Sampling focal coordinates for attempt %d...\n", attempt))
    coords <- .synthetic_helpers$sample_focal_coords(
      raster_bundle$rast,
      .synthetic_helpers$candidate_population_count(base_config),
      seed = base_config$seed + as.integer(config$dim) + 1000L * as.integer(attempt)
    )
    saveRDS(coords, coords_path)
  }

  if (file.exists(cost_path))
  {
    cost_bundle <- readRDS(cost_path)
    cost.mat <- cost_bundle$cost.mat
  }
  else
  {
    cat(sprintf("Computing cost matrix for attempt %d...\n", attempt))
    surface_candidate <- terradish::conductance_surface(
      raster_bundle$covariates,
      coords,
      directions = base_config$directions
    )
    cost.mat <- .synthetic_helpers$true_cost_matrix(surface_candidate, base_config)
    saveRDS(list(cost.mat = cost.mat), cost_path)
  }

  if (file.exists(fst_path))
  {
    fst <- readRDS(fst_path)
  }
  else
  {
    cat(sprintf("Simulating FST for attempt %d...\n", attempt))
    fst <- .synthetic_helpers$simulate_fst(coords, cost.mat, base_config)
    saveRDS(fst, fst_path)
  }
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
  attempt_result <- list(
    attempt = as.integer(attempt),
    coords = coords,
    cost.mat = cost.mat,
    fst = fst,
    n_surviving = nrow(coords)
  )
  saveRDS(attempt_result, attempt_path)
  attempt_result
}

.load_surface_case <- function(paths, config, base_config)
{
  if (file.exists(paths$surface_case))
    return(.unpack_surface_case(readRDS(paths$surface_case)))

  raster_bundle <- .load_surface_raster(paths, config, base_config)

  selected_bundle <- NULL
  if (file.exists(paths$surface_selected))
    selected_bundle <- readRDS(paths$surface_selected)

  if (is.null(selected_bundle))
  {
    for (attempt in seq_len(base_config$max_sim_attempts))
    {
      attempt_result <- .load_surface_attempt(
        paths = paths,
        output_dir = config$output_dir,
        attempt = attempt,
        raster_bundle = raster_bundle,
        config = config,
        base_config = base_config
      )
      if (attempt_result$n_surviving >= base_config$n_pops)
      {
        selected <- sort(sample(seq_len(attempt_result$n_surviving),
                                size = base_config$n_pops,
                                replace = FALSE))
        selected_bundle <- list(
          attempt = attempt_result$attempt,
          coords = attempt_result$coords[selected, , drop = FALSE],
          fst = attempt_result$fst[selected, selected, drop = FALSE],
          selected = selected
        )
        saveRDS(selected_bundle, paths$surface_selected)
        break
      }
    }
  }

  if (is.null(selected_bundle))
    stop("Genetic simulation returned fewer than the requested number of focal populations")
  if (nrow(selected_bundle$coords) != nrow(selected_bundle$fst))
    stop("Simulated focal coordinates and FST matrix are misaligned")

  cat("Building final conductance surface from selected focal populations...\n")
  surface <- terradish::conductance_surface(
    raster_bundle$covariates,
    selected_bundle$coords,
    directions = base_config$directions
  )
  surface_case <- list(
    rast = raster_bundle$rast,
    covariates = raster_bundle$covariates,
    coords = selected_bundle$coords,
    fst = selected_bundle$fst,
    surface = surface,
    selected_attempt = selected_bundle$attempt
  )
  saveRDS(.pack_surface_case(surface_case), paths$surface_case)
  surface_case
}

.load_exact_fit <- function(paths, surface_case, base_config)
{
  if (file.exists(paths$exact_fit))
    return(readRDS(paths$exact_fit))

  cat("Running exact fit baseline...\n")
  fst <- surface_case$fst
  exact_fit <- .elapsed_run(suppressWarnings(
    terradish::radish(
      fst ~ cont1 + cont2,
      data = surface_case$surface,
      conductance_model = terradish::loglinear_conductance,
      measurement_model = terradish::generalized_wishart,
      nu = base_config$n_loci,
      control = terradish::NewtonRaphsonControl(maxit = 6, verbose = FALSE),
      solver = "direct"
    )
  ))
  saveRDS(exact_fit, paths$exact_fit)
  exact_fit
}

.run_exact_grid <- function(theta, surface_case, base_config)
{
  fst <- surface_case$fst
  .elapsed_run(terradish::radish_grid(
    theta = theta,
    formula = fst ~ cont1 + cont2,
    data = surface_case$surface,
    conductance_model = terradish::loglinear_conductance,
    measurement_model = terradish::generalized_wishart,
    nu = base_config$n_loci,
    cores = 1,
    approximation = "none"
  ))
}

.run_coarse_grid <- function(theta, surface_case, base_config, factor)
{
  fst <- surface_case$fst
  .elapsed_run(terradish::radish_grid(
    theta = theta,
    formula = fst ~ cont1 + cont2,
    data = surface_case$surface,
    conductance_model = terradish::loglinear_conductance,
    measurement_model = terradish::generalized_wishart,
    nu = base_config$n_loci,
    cores = 1,
    approximation = "coarse_raster",
    approximation_control = list(factor = as.integer(factor))
  ))
}

.run_exact_fit_from_start <- function(surface_case, base_config, theta_start)
{
  fst <- surface_case$fst
  .elapsed_run(suppressWarnings(
    terradish::radish(
      fst ~ cont1 + cont2,
      data = surface_case$surface,
      conductance_model = terradish::loglinear_conductance,
      measurement_model = terradish::generalized_wishart,
      nu = base_config$n_loci,
      theta = theta_start,
      control = terradish::NewtonRaphsonControl(maxit = 6, verbose = FALSE),
      solver = "direct"
    )
  ))
}

main <- function(config = CONFIG)
{
  if (!requireNamespace("terradish", quietly = TRUE))
    stop("Package `terradish` must be installed.")
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package `terra` must be installed.")
  if (!requireNamespace("multiScaleR", quietly = TRUE))
    stop("Package `multiScaleR` must be installed.")
  if (!requireNamespace("PopGenReport", quietly = TRUE))
    stop("Package `PopGenReport` must be installed.")

  config$output_dir <- .ensure_output_dir(config$output_dir)
  path_map <- .paths(config$output_dir)
  base_config <- .merge_screening_config(config)

  surface_case <- .load_surface_case(path_map, config, base_config)

  center <- switch(
    config$center_mode,
    true_theta = base_config$true_theta,
    exact_fit = .load_exact_fit(path_map, surface_case, base_config)$value$mle$theta,
    stop("Unknown `center_mode`: ", config$center_mode)
  )
  theta <- .build_theta_grid(center,
                             half_width = config$theta_half_width,
                             grid_points = as.integer(config$grid_points))

  grid_results <- .read_or_empty(
    path_map$grid_raw,
    prototype = data.frame(
      mode = character(),
      elapsed_sec = numeric(),
      max_logLik = numeric(),
      min_logLik = numeric(),
      spearman = numeric(),
      max_abs_logLik_diff = numeric(),
      mean_abs_logLik_diff = numeric(),
      top_match = logical(),
      approx_type = character(),
      approx_used = logical(),
      full_vertices = numeric(),
      coarse_vertices = numeric(),
      stringsAsFactors = FALSE
    )
  )

  if (!"exact" %in% grid_results$mode)
  {
    cat("Running exact grid...\n")
    exact_run <- .run_exact_grid(theta, surface_case, base_config)
    exact_row <- .grid_result_row("exact", exact_run, surface_case, baseline = NULL)
    grid_results <- rbind(grid_results, exact_row)
    saveRDS(exact_run, path_map$exact_grid)
    .write_csv(.summarize_grid_results(grid_results), path_map$grid_summary)
    .write_csv(grid_results, path_map$grid_raw)
  }

  if (file.exists(path_map$exact_grid))
    exact_grid_result <- readRDS(path_map$exact_grid)$value
  else
    exact_grid_result <- .run_exact_grid(theta, surface_case, base_config)$value

  for (factor in as.integer(config$factors))
  {
    mode_name <- paste0("coarse_factor_", factor)
    if (mode_name %in% grid_results$mode)
      next

    cat(sprintf("Running %s...\n", mode_name))
    run <- .run_coarse_grid(theta, surface_case, base_config, factor)
    row <- .grid_result_row(mode_name, run, surface_case, baseline = exact_grid_result)
    grid_results <- rbind(grid_results, row)
    saveRDS(run, .grid_rds_path(config$output_dir, mode_name))
    .write_csv(grid_results, path_map$grid_raw)
    .write_csv(.summarize_grid_results(grid_results), path_map$grid_summary)
  }

  grid_summary <- .summarize_grid_results(grid_results)

  fit_summary <- NULL
  if (isTRUE(config$benchmark_starts))
  {
    exact_fit <- .load_exact_fit(path_map, surface_case, base_config)
    fit_results <- .read_or_empty(
      path_map$start_raw,
      prototype = data.frame(
        start = character(),
        elapsed_sec = numeric(),
        logLik = numeric(),
        max_abs_theta_diff = numeric(),
        stringsAsFactors = FALSE
      )
    )

    if (!"default" %in% fit_results$start)
    {
      fit_results <- rbind(fit_results, .fit_result_row("default", exact_fit, exact_fit$value))
      .write_csv(fit_results, path_map$start_raw)
      .write_csv(.summarize_fit_results(fit_results), path_map$start_summary)
    }

    theta_df <- theta
    for (factor in as.integer(config$factors))
    {
      label <- paste0("factor_", factor)
      if (label %in% fit_results$start)
        next

      grid_row <- grid_results[grid_results$mode == paste0("coarse_factor_", factor), , drop = FALSE]
      if (!nrow(grid_row))
        next

      coarse_path <- .grid_rds_path(config$output_dir, paste0("coarse_factor_", factor))
      coarse_grid <- if (file.exists(coarse_path))
        readRDS(coarse_path)$value
      else
        .run_coarse_grid(theta, surface_case, base_config, factor)$value
      best_idx <- which.max(coarse_grid$loglik)
      theta_start <- as.numeric(theta_df[best_idx, , drop = TRUE])
      names(theta_start) <- colnames(theta_df)

      cat(sprintf("Running exact fit from %s grid winner...\n", label))
      run <- .run_exact_fit_from_start(surface_case, base_config, theta_start)
      fit_results <- rbind(fit_results, .fit_result_row(label, run, exact_fit$value))
      .write_csv(fit_results, path_map$start_raw)
      .write_csv(.summarize_fit_results(fit_results), path_map$start_summary)
    }

    fit_summary <- .summarize_fit_results(fit_results)
  }

  cat("\nGrid summary:\n")
  print(grid_summary, row.names = FALSE)
  if (!is.null(fit_summary))
  {
    cat("\nFit-start summary:\n")
    print(fit_summary, row.names = FALSE)
  }

  invisible(list(
    grid = grid_summary,
    fit_starts = fit_summary,
    output_dir = config$output_dir
  ))
}

if (sys.nframe() == 0L)
  main()
