#!/usr/bin/env Rscript

# Focused curvature benchmark for the large-landscape workflow.
#
# This script compares curvature = "exact" and curvature = "gauss_newton" on a
# synthetic heterogeneous raster with a generalized Wishart distance response.
# The default problem is intentionally modest so the script is reviewable and
# can run on a laptop. Increase CONFIG$side and set CONFIG$solver = "amg" to
# repeat the same comparison in the million-cell large-N regime.
#
# Usage:
#   source("inst/benchmarks/benchmark-curvature-large-landscape.R")
# or
#   Rscript inst/benchmarks/benchmark-curvature-large-landscape.R

CONFIG <- list(
  seed = 11L,
  side = 96L,
  n_focal = 24L,
  directions = 8L,
  true_theta = c(ridge = 0.6, basin = -0.45),
  tau = 1,
  sigma = 0.05,
  nu = 1000L,
  solver = "auto",
  solver_control = NULL,
  optimizer = "bfgs",
  maxit = 50L,
  output_csv = NULL
)

.require_packages <- function()
{
  needed <- c("terra", "terradish")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

.elapsed <- function(expr)
{
  gc()
  started <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = unname(proc.time()[["elapsed"]] - started))
}

.make_layer <- function(values, side)
{
  r <- terra::rast(nrows = side, ncols = side, xmin = 0, xmax = side,
                   ymin = 0, ymax = side)
  terra::values(r) <- as.vector(values)
  r
}

.make_synthetic_case <- function(config)
{
  side <- as.integer(config$side)
  set.seed(as.integer(config$seed))

  x <- seq(0, 1, length.out = side)
  y <- seq(0, 1, length.out = side)
  gx <- matrix(rep(x, each = side), nrow = side)
  gy <- matrix(rep(y, times = side), nrow = side)

  ridge <- sin(3 * pi * gx) + 0.5 * cos(2 * pi * gy) +
    0.20 * matrix(stats::rnorm(side * side), side, side)
  basin <- exp(-((gx - 0.62)^2 + (gy - 0.38)^2) / 0.035) -
    0.35 * exp(-((gx - 0.24)^2 + (gy - 0.76)^2) / 0.020) +
    0.15 * matrix(stats::rnorm(side * side), side, side)

  covariates <- c(.make_layer(ridge, side), .make_layer(basin, side))
  names(covariates) <- names(config$true_theta)
  covariates <- terradish::scale_covariates(covariates)

  focal_cells <- sort(sample.int(terra::ncell(covariates[[1]]),
                                 as.integer(config$n_focal)))
  coords <- terra::vect(
    terra::xyFromCell(covariates[[1]], focal_cells),
    type = "points",
    crs = terra::crs(covariates)
  )

  surface <- terradish::conductance_surface(
    covariates,
    coords,
    directions = as.integer(config$directions),
    saveStack = TRUE
  )

  theta <- matrix(unname(config$true_theta), nrow = 1L)
  colnames(theta) <- names(config$true_theta)
  sim <- terradish::simulate_covariance_response(
    theta = theta,
    formula = stats::reformulate(names(config$true_theta)),
    data = surface,
    conductance_model = terradish::loglinear_conductance,
    tau = config$tau,
    sigma = config$sigma,
    nu = as.integer(config$nu),
    seed = as.integer(config$seed) + 1000L
  )

  list(
    surface = surface,
    response = terradish::dist_from_cov(sim$covariance),
    formula = stats::as.formula(
      paste("S_dist ~", paste(names(config$true_theta), collapse = " + "))
    ),
    theta_start = stats::setNames(rep(0, length(config$true_theta)),
                                  names(config$true_theta)),
    n_vertices = nrow(surface$x),
    n_focal = nrow(sim$covariance)
  )
}

.fit_once <- function(case, curvature, config)
{
  S_dist <- case$response
  fit_formula <- case$formula
  environment(fit_formula) <- environment()

  control <- terradish::NewtonRaphsonControl(
    maxit = as.integer(config$maxit),
    verbose = FALSE
  )

  fit <- terradish::terradish(
    fit_formula,
    data = case$surface,
    conductance_model = terradish::loglinear_conductance,
    measurement_model = terradish::generalized_wishart,
    nu = as.integer(config$nu),
    theta = case$theta_start,
    optimizer = config$optimizer,
    solver = config$solver,
    solver_control = config$solver_control,
    control = control,
    curvature = curvature
  )

  summary_fit <- summary(fit)
  hessian_eigen <- eigen(fit$fit$hessian, symmetric = TRUE,
                         only.values = TRUE)$values
  solver_info <- fit$fit$solver_info

  data.frame(
    curvature = curvature,
    elapsed_sec = NA_real_,
    logLik = fit$loglik,
    aic = fit$aic,
    gradient_norm = sqrt(sum(fit$mle$gradient^2)),
    min_hessian_eigen = min(hessian_eigen),
    max_hessian_eigen = max(hessian_eigen),
    hessian_condition = max(abs(hessian_eigen)) / max(min(abs(hessian_eigen)),
                                                       .Machine$double.eps),
    solver = if (is.null(solver_info$type)) NA_character_ else solver_info$type,
    auto_reason = if (is.null(solver_info$auto_reason)) "" else solver_info$auto_reason,
    function_calls = unname(fit$cost[["function_calls"]]),
    optimizer_steps = unname(fit$cost[["newton_steps"]]),
    coef_ridge = unname(stats::coef(fit)[["ridge"]]),
    coef_basin = unname(stats::coef(fit)[["basin"]]),
    se_ridge = unname(summary_fit$ztable["ridge", "Std. Error"]),
    se_basin = unname(summary_fit$ztable["basin", "Std. Error"]),
    boundary = isTRUE(fit$fit$boundary),
    stringsAsFactors = FALSE
  )
}

.summarize_pair <- function(results)
{
  exact <- results[results$curvature == "exact", , drop = FALSE]
  gn <- results[results$curvature == "gauss_newton", , drop = FALSE]
  data.frame(
    max_abs_coef_diff = max(abs(c(exact$coef_ridge - gn$coef_ridge,
                                  exact$coef_basin - gn$coef_basin))),
    abs_logLik_diff = abs(exact$logLik - gn$logLik),
    ridge_se_ratio_gn_over_exact = gn$se_ridge / exact$se_ridge,
    basin_se_ratio_gn_over_exact = gn$se_basin / exact$se_basin,
    elapsed_ratio_gn_over_exact = gn$elapsed_sec / exact$elapsed_sec,
    exact_min_hessian_eigen = exact$min_hessian_eigen,
    gn_min_hessian_eigen = gn$min_hessian_eigen,
    stringsAsFactors = FALSE
  )
}

.print_recommendation <- function(results, comparison)
{
  cat("\nRecommendation:\n")
  cat("- Use curvature = \"exact\" for routine final fits when the Hessian is\n")
  cat("  positive definite and runtime is acceptable; it is the observed-curvature\n")
  cat("  default and is the most literal estimate of local likelihood shape.\n")
  cat("- Use curvature = \"gauss_newton\" in the large-landscape workflow when you\n")
  cat("  want positive-semidefinite, information-based standard errors, or when the\n")
  cat("  exact Hessian is indefinite or unstable near the fit.\n")
  cat("- Re-run both settings when the Gauss-Newton and exact standard errors differ\n")
  cat("  materially; a large gap means the residual-weighted second-derivative terms\n")
  cat("  are carrying information, often because the fit is not close to the model\n")
  cat("  mean or the model is misspecified.\n")
  cat(sprintf("- In this run, max |coefficient difference| was %.3g and the\n",
              comparison$max_abs_coef_diff))
  cat(sprintf("  log-likelihood difference was %.3g, so the curvature choice did\n",
              comparison$abs_logLik_diff))
  cat("  not move the fitted surface.\n")
}

run_curvature_large_landscape_benchmark <- function(config = CONFIG)
{
  .require_packages()
  case <- .make_synthetic_case(config)

  cat(sprintf("Synthetic case: %d x %d raster, %d retained vertices, %d focal sites\n",
              as.integer(config$side), as.integer(config$side),
              case$n_vertices, case$n_focal))
  cat(sprintf("Solver: %s; optimizer: %s; nu: %d\n",
              config$solver, config$optimizer, as.integer(config$nu)))

  rows <- lapply(c("exact", "gauss_newton"), function(curvature) {
    cat(sprintf("  fitting curvature = %s\n", curvature))
    timed <- .elapsed(.fit_once(case, curvature, config))
    timed$value$elapsed_sec <- timed$elapsed
    timed$value
  })

  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  comparison <- .summarize_pair(results)

  cat("\nRaw results:\n")
  print(results, row.names = FALSE)
  cat("\nPairwise comparison:\n")
  print(comparison, row.names = FALSE)
  .print_recommendation(results, comparison)

  if (!is.null(config$output_csv))
    utils::write.csv(results, config$output_csv, row.names = FALSE)

  invisible(list(results = results, comparison = comparison, config = config))
}

main <- function(config = CONFIG)
  run_curvature_large_landscape_benchmark(config)

if (sys.nframe() == 0L)
  main(CONFIG)
