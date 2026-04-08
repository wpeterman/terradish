#!/usr/bin/env Rscript

# Focused speed benchmark for terradish solver/optimizer settings.
#
# This script compares a small, reviewable matrix of settings on:
# - the built-in melip example
# - a synthetic 250 x 250 raster graph
# - a synthetic 750 x 750 raster graph
#
# The goal is not exhaustive tuning. Instead, it provides a reproducible
# benchmark for the settings that have been most relevant in development:
# direct sparse solves, AMG, and coarse-raster refinement.
#
# Baseline observations from the April 2026 exploratory run on the development
# machine:
# - `direct` was fastest on all tested problems.
# - Explicit `amg` reached the same objective values, but was slower on melip,
#   on a synthetic 250 x 250 grid, and on a synthetic 750 x 750 grid.
# - `coarse_raster` was a useful option to keep benchmarking, but it did not
#   beat the best exact direct fit once exact refinement was included.
# - For the larger synthetic problem tested, `supernodal_ll` was the direct
#   factorization selected and remained competitive.
#
# Treat those as a starting point, not a guarantee: rerun this script after
# solver, optimizer, or compiled-code changes.
#
# Usage:
#   source("inst/benchmarks/benchmark-speed-settings.R")
# or
#   Rscript inst/benchmarks/benchmark-speed-settings.R

CONFIG <- list(
  seed = 1L,
  output_csv = NULL,
  synthetic_dims = c(250L, 750L),
  synthetic_n_covariates = 5L,
  synthetic_n_focal = 20L,
  synthetic_wishart_df = 60L,
  settings = list(
    melip = list(
      list(label = "direct_newton", optimizer = "newton", solver = "direct",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = TRUE),
      list(label = "direct_bfgs", optimizer = "bfgs", solver = "direct",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = TRUE),
      list(label = "amg_newton", optimizer = "newton", solver = "amg",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = TRUE),
      list(label = "coarse_auto", optimizer = "auto", solver = "direct",
           approximation = "coarse_raster",
           approximation_control = list(factor = 2L, exact_refine = TRUE),
           measurement_model = "mlpe", nonnegative = TRUE)
    ),
    synthetic = list(
      list(label = "direct_newton", optimizer = "newton", solver = "direct",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = FALSE),
      list(label = "direct_bfgs", optimizer = "bfgs", solver = "direct",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = FALSE),
      list(label = "amg_bfgs", optimizer = "bfgs", solver = "amg",
           approximation = "none", approximation_control = NULL,
           measurement_model = "mlpe", nonnegative = FALSE),
      list(label = "coarse_auto", optimizer = "auto", solver = "direct",
           approximation = "coarse_raster",
           approximation_control = list(factor = 4L, exact_refine = TRUE),
           measurement_model = "mlpe", nonnegative = FALSE)
    )
  )
)

.require_packages <- function()
{
  needed <- c("terra", "terradish")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

.elapsed_run <- function(expr)
{
  gc()
  started <- Sys.time()
  result <- force(expr)
  elapsed <- unname((proc.time())[["elapsed"]])
  list(
    value = result,
    started = started,
    elapsed = elapsed
  )
}

.timed_fit <- function(case, setting)
{
  fit_formula <- case$formula
  had_response <- exists(case$response_name, envir = .GlobalEnv, inherits = FALSE)
  if (had_response)
    old_response <- get(case$response_name, envir = .GlobalEnv, inherits = FALSE)
  assign(case$response_name, case$response, envir = .GlobalEnv)
  on.exit({
    if (had_response)
      assign(case$response_name, old_response, envir = .GlobalEnv)
    else if (exists(case$response_name, envir = .GlobalEnv, inherits = FALSE))
      rm(list = case$response_name, envir = .GlobalEnv)
  }, add = TRUE)
  environment(fit_formula) <- .GlobalEnv

  control <- terradish::NewtonRaphsonControl(
    verbose = FALSE,
    ctol = 1e-6,
    ftol = 1e-6
  )

  started <- Sys.time()
  timer <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    suppressWarnings(
      terradish::terradish(
        fit_formula,
        data = case$surface,
        conductance_model = terradish::loglinear_conductance,
        measurement_model = get(setting$measurement_model, envir = asNamespace("terradish")),
        theta = case$theta_start,
        nonnegative = setting$nonnegative,
        optimizer = setting$optimizer,
        solver = setting$solver,
        approximation = setting$approximation,
        approximation_control = setting$approximation_control,
        control = control
      )
    ),
    error = function(e) e
  )
  elapsed <- unname((proc.time()[["elapsed"]] - timer))
  ended <- Sys.time()

  if (inherits(fit, "error"))
  {
    return(data.frame(
      dataset = case$name,
      setting = setting$label,
      optimizer = setting$optimizer,
      solver = setting$solver,
      approximation = setting$approximation,
      elapsed = elapsed,
      loglik = NA_real_,
      steps = NA_integer_,
      calls = NA_integer_,
      solver_used = NA_character_,
      factorization = NA_character_,
      auto_reason = NA_character_,
      approx_used = NA,
      coarse_vertices = NA_integer_,
      full_vertices = nrow(case$surface$x),
      status = paste("ERROR:", conditionMessage(fit)),
      started = as.character(started),
      ended = as.character(ended),
      stringsAsFactors = FALSE
    ))
  }

  sinfo <- fit$fit$solver_info
  approx <- fit$approximation

  data.frame(
    dataset = case$name,
    setting = setting$label,
    optimizer = setting$optimizer,
    solver = setting$solver,
    approximation = setting$approximation,
    elapsed = elapsed,
    loglik = fit$loglik,
    steps = .named_cost_value(fit$cost, "newton_steps"),
    calls = .named_cost_value(fit$cost, "function_calls"),
    solver_used = if (is.null(sinfo$type)) NA_character_ else sinfo$type,
    factorization = if (is.null(sinfo$factorization)) NA_character_ else sinfo$factorization,
    auto_reason = if (is.null(sinfo$auto_reason)) NA_character_ else sinfo$auto_reason,
    approx_used = isTRUE(approx$used),
    coarse_vertices = if (is.null(approx$coarse_vertices)) NA_integer_ else approx$coarse_vertices,
    full_vertices = nrow(case$surface$x),
    status = "OK",
    started = as.character(started),
    ended = as.character(ended),
    stringsAsFactors = FALSE
  )
}

.print_progress <- function(i, n, dataset, label)
{
  cat(sprintf("[%02d/%02d] %s :: %s\n", i, n, dataset, label))
  flush.console()
}

.named_cost_value <- function(cost, key)
{
  nms <- names(cost)
  if (is.null(nms) || !(key %in% nms))
    return(NA_integer_)
  unname(as.integer(cost[[key]]))
}

.summarize_results <- function(results)
{
  split_results <- split(results, results$dataset)
  do.call(rbind, lapply(split_results, function(df)
  {
    ordered <- df[order(df$elapsed), , drop = FALSE]
    ordered$speedup_vs_fastest <- ordered$elapsed[[1]] / ordered$elapsed
    ordered
  }))
}

.make_melip_case <- function()
{
  data("melip", package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- terradish::scale_covariates(covariates)

  list(
    name = "melip",
    response_name = "melip.Fst",
    response = melip.Fst,
    formula = melip.Fst ~ forestcover + altitude,
    surface = terradish::conductance_surface(covariates, melip.coords,
                                             directions = 8L, saveStack = TRUE),
    theta_start = NULL
  )
}

.make_synthetic_case <- function(dim, config)
{
  k <- as.integer(config$synthetic_n_covariates)
  n_focal <- as.integer(config$synthetic_n_focal)
  seed <- as.integer(config$seed + dim)

  set.seed(seed)
  xseq <- seq(0, 1, length.out = dim)
  yseq <- seq(0, 1, length.out = dim)
  grid_x <- matrix(rep(xseq, each = dim), nrow = dim, ncol = dim)
  grid_y <- matrix(rep(yseq, times = dim), nrow = dim, ncol = dim)

  layers <- vector("list", k)
  for (j in seq_len(k))
  {
    signal <- sin(pi * j * grid_x) +
      cos(pi * (j + 1) * grid_y) +
      0.35 * sin(2 * pi * (grid_x + grid_y) / (j + 1))
    noise <- matrix(rnorm(dim * dim, sd = 0.10), nrow = dim, ncol = dim)

    r <- terra::rast(nrows = dim, ncols = dim, xmin = 0, xmax = 1,
                     ymin = 0, ymax = 1)
    terra::values(r) <- as.vector(signal + noise)
    layers[[j]] <- r
  }

  covariates <- do.call(c, layers)
  names(covariates) <- paste0("var", seq_len(k))
  covariates <- terradish::scale_covariates(covariates)

  focal_cells <- sort(sample.int(terra::ncell(covariates[[1]]), n_focal))
  coords <- terra::vect(
    terra::xyFromCell(covariates[[1]], focal_cells),
    type = "points",
    crs = terra::crs(covariates)
  )

  surface <- terradish::conductance_surface(
    covariates,
    coords,
    directions = 4L,
    saveStack = TRUE
  )

  theta_true <- matrix(c(1.2, -1.0, 0.8, -0.7, 0.5), nrow = 1)
  colnames(theta_true) <- paste0("var", seq_len(k))

  had_loglinear <- exists("loglinear_conductance", envir = .GlobalEnv, inherits = FALSE)
  if (!had_loglinear)
    assign("loglinear_conductance", terradish::loglinear_conductance, envir = .GlobalEnv)
  on.exit({
    if (!had_loglinear && exists("loglinear_conductance", envir = .GlobalEnv, inherits = FALSE))
      rm("loglinear_conductance", envir = .GlobalEnv)
  }, add = TRUE)

  covariance <- terradish::terradish_distance(
    theta = theta_true,
    formula = stats::reformulate(paste0("var", seq_len(k))),
    data = surface,
    conductance_model = terradish::loglinear_conductance,
    conductance = TRUE,
    covariance = TRUE
  )$covariance[, , 1]

  set.seed(seed + 100L)
  response <- terradish::dist_from_cov(
    solve(stats::rWishart(1, as.integer(config$synthetic_wishart_df), covariance)[,,1])
  )

  list(
    name = sprintf("synthetic_%dx%d", dim, dim),
    response_name = "S",
    response = response,
    formula = stats::as.formula(
      paste("S ~", paste(paste0("var", seq_len(k)), collapse = " + "))
    ),
    surface = surface,
    theta_start = rep(0, k)
  )
}

run_speed_settings_benchmark <- function(config = CONFIG)
{
  .require_packages()

  cases <- c(
    list(.make_melip_case()),
    lapply(config$synthetic_dims, function(dim) .make_synthetic_case(as.integer(dim), config))
  )

  settings <- c(
    lapply(config$settings$melip, function(x) list(case = "melip", setting = x)),
    unlist(lapply(sprintf("synthetic_%dx%d", config$synthetic_dims, config$synthetic_dims),
                  function(name) {
                    lapply(config$settings$synthetic, function(x) list(case = name, setting = x))
                  }), recursive = FALSE)
  )

  case_map <- stats::setNames(cases, vapply(cases, `[[`, character(1), "name"))
  out <- vector("list", length(settings))

  for (i in seq_along(settings))
  {
    spec <- settings[[i]]
    case <- case_map[[spec$case]]
    setting <- spec$setting
    .print_progress(i, length(settings), case$name, setting$label)
    out[[i]] <- .timed_fit(case, setting)
    print(out[[i]][, c("dataset", "setting", "elapsed", "loglik", "steps",
                       "calls", "solver_used", "factorization",
                       "approx_used", "status")])
    flush.console()
  }

  results <- do.call(rbind, out)
  ranked <- .summarize_results(results)

  if (!is.null(config$output_csv))
  {
    utils::write.csv(results, config$output_csv, row.names = FALSE)
  }

  list(results = results, ranked = ranked)
}

main <- function(config = CONFIG)
{
  benchmark <- run_speed_settings_benchmark(config)

  cat("\n=== BENCHMARK RESULTS ===\n")
  print(benchmark$results)

  cat("\n=== RANKED BY ELAPSED WITHIN DATASET ===\n")
  print(benchmark$ranked[, c("dataset", "setting", "elapsed", "loglik", "steps",
                             "calls", "solver_used", "factorization",
                             "approx_used", "speedup_vs_fastest", "status")])

  invisible(benchmark)
}

if (sys.nframe() == 0L)
  main(CONFIG)
