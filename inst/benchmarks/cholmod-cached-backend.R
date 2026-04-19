#!/usr/bin/env Rscript

# Focused benchmark for terradish direct CHOLMOD backends.
#
# This compares:
# - "matrix": Matrix::Cholesky/update/solve through the standard R interface
# - "cholmod_cpp": Matrix factorization/update plus compiled CHOLMOD solve
# - "cholmod_cpp_cached": compiled CHOLMOD factorization/update/solve held in
#   an external pointer across nearby evaluations
#
# The key timing is the second "nearby" evaluation because that is the pattern
# used during optimization and line search. First evaluations include symbolic
# analysis, while second evaluations can reuse compatible solver state.
#
# Usage from the package root:
#   Rscript inst/benchmarks/cholmod-cached-backend.R
#
# Optional environment variables:
#   TERRADISH_BENCH_GRIDS="50,100,200,300"
#   TERRADISH_BENCH_REPS=3
#   TERRADISH_BENCH_OUTPUT="cholmod-backend-results.csv"
#   TERRADISH_BENCH_LOAD_ALL=false

CONFIG <- list(
  grids = c(50L, 100L, 200L, 300L),
  n_sites = 25L,
  directions = 8L,
  reps = 3L,
  backends = c("matrix", "cholmod_cpp", "cholmod_cpp_cached"),
  factorization = "auto",
  output_csv = Sys.getenv("TERRADISH_BENCH_OUTPUT", unset = NA_character_),
  load_all = !identical(tolower(Sys.getenv("TERRADISH_BENCH_LOAD_ALL", "true")), "false")
)

parse_int_env <- function(name, default)
{
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value))
    return(default)
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  out <- out[!is.na(out)]
  if (!length(out))
    return(default)
  out
}

load_terradish_for_benchmark <- function()
{
  if (isTRUE(CONFIG$load_all) &&
      file.exists("DESCRIPTION") &&
      requireNamespace("pkgload", quietly = TRUE))
  {
    pkgload::load_all(".", export_all = TRUE, helpers = FALSE)
    return(invisible(TRUE))
  }

  library(terradish)
  invisible(TRUE)
}

make_melip_surface <- function()
{
  data(melip, package = "terradish", envir = environment())
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- scale_covariates(covariates)

  surface <- conductance_surface(covariates, melip.coords, directions = CONFIG$directions)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)

  list(
    label = "melip",
    surface = surface,
    C1 = model(c(altitude = -0.3, forestcover = 0.3))$conductance,
    C2 = model(c(altitude = -0.28, forestcover = 0.32))$conductance
  )
}

make_synthetic_surface <- function(n, n_sites = CONFIG$n_sites)
{
  r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- scale(xy[, 1] / n + 0.5 * xy[, 2] / n)[, 1]
  names(r) <- "gradient"

  side <- ceiling(sqrt(n_sites))
  rows <- unique(pmin(pmax(round(seq(2, n - 1, length.out = side)), 1L), n))
  cols <- unique(pmin(pmax(round(seq(2, n - 1, length.out = side)), 1L), n))
  grid <- expand.grid(row = rows, col = cols)
  cells <- terra::cellFromRowCol(r, grid$row, grid$col)[seq_len(n_sites)]
  coords <- terra::vect(
    as.data.frame(terra::xyFromCell(r, cells)),
    geom = c("x", "y"),
    crs = terra::crs(r)
  )

  surface <- conductance_surface(r, coords, directions = CONFIG$directions)
  x <- as.numeric(surface$x[, 1])

  list(
    label = paste0(n, "x", n),
    surface = surface,
    C1 = exp(0.15 * x),
    C2 = exp(0.18 * x)
  )
}

solver_reuse_state <- function(state, backend)
{
  if (identical(backend, "cholmod_cpp_cached"))
  {
    return(list(
      type = "direct",
      handle = state$handle,
      signature = state$signature
    ))
  }

  list(
    type = "direct",
    factor = state$factor,
    signature = state$signature
  )
}

benchmark_one_backend <- function(case, backend, rep_id, ref_update = NULL)
{
  surface <- case$surface
  rhs <- terradish:::.graph_rhs(surface, length(case$C1))
  control <- list(factorization = CONFIG$factorization, solve_backend = backend)

  gc(FALSE)
  first_start <- proc.time()[["elapsed"]]
  state1 <- terradish:::.terradish_solver_setup(
    surface,
    case$C1,
    solver = "direct",
    solver_control = control
  )
  solve1 <- terradish:::.terradish_solver_solve(state1, rhs)
  first_total <- proc.time()[["elapsed"]] - first_start

  gc(FALSE)
  update_start <- proc.time()[["elapsed"]]
  state2 <- terradish:::.terradish_solver_setup(
    surface,
    case$C2,
    solver = "direct",
    solver_control = control,
    solver_reuse_state = solver_reuse_state(state1, backend)
  )
  solve2 <- terradish:::.terradish_solver_solve(state2, rhs)
  update_total <- proc.time()[["elapsed"]] - update_start

  update_solution <- as.matrix(solve2$solution)
  max_abs_diff <- if (is.null(ref_update) && identical(backend, "matrix"))
    0
  else if (is.null(ref_update))
    NA_real_
  else
    max(abs(update_solution - ref_update))

  list(
    row = data.frame(
      dataset = case$label,
      backend = backend,
      replicate = rep_id,
      n_vertices = nrow(surface$x),
      n_rhs = ncol(rhs),
      factorization = solve2$info$factorization,
      first_total = first_total,
      first_setup = unname(solve1$info$setup_time),
      first_solve = unname(solve1$info$solve_time),
      update_total = update_total,
      update_setup = unname(solve2$info$setup_time),
      update_solve = unname(solve2$info$solve_time),
      reused = isTRUE(solve2$info$reused_factor_template),
      max_abs_diff = max_abs_diff,
      stringsAsFactors = FALSE
    ),
    update_solution = update_solution
  )
}

benchmark_case <- function(case)
{
  rows <- list()
  row_id <- 0L

  for (rep_id in seq_len(CONFIG$reps))
  {
    ref_update <- NULL
    for (backend in CONFIG$backends)
    {
      out <- benchmark_one_backend(case, backend, rep_id, ref_update = ref_update)
      if (identical(backend, "matrix"))
        ref_update <- out$update_solution

      row_id <- row_id + 1L
      rows[[row_id]] <- out$row
    }
  }

  do.call(rbind, rows)
}

summarize_results <- function(results)
{
  summary <- stats::aggregate(
    cbind(first_total, first_setup, first_solve,
          update_total, update_setup, update_solve, max_abs_diff) ~
      dataset + n_vertices + n_rhs + factorization + backend,
    results,
    stats::median,
    na.rm = TRUE
  )

  summary[order(summary$n_vertices, summary$backend), , drop = FALSE]
}

main <- function()
{
  CONFIG$grids <<- parse_int_env("TERRADISH_BENCH_GRIDS", CONFIG$grids)
  CONFIG$reps <<- parse_int_env("TERRADISH_BENCH_REPS", CONFIG$reps)[1]

  load_terradish_for_benchmark()

  cases <- list(make_melip_surface())
  for (n in CONFIG$grids)
  {
    message("Building synthetic ", n, "x", n)
    cases[[length(cases) + 1L]] <- make_synthetic_surface(n)
  }

  results <- list()
  for (case in cases)
  {
    message("Benchmarking ", case$label)
    results[[length(results) + 1L]] <- benchmark_case(case)
  }

  results <- do.call(rbind, results)
  summary <- summarize_results(results)

  print(summary, row.names = FALSE, digits = 4)

  if (!is.na(CONFIG$output_csv) && nzchar(CONFIG$output_csv))
  {
    utils::write.csv(results, CONFIG$output_csv, row.names = FALSE)
    message("Wrote raw benchmark results to ", CONFIG$output_csv)
  }

  invisible(list(raw = results, summary = summary))
}

if (identical(environment(), globalenv()))
  main()
