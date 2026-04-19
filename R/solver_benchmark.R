#' Benchmark linear solver settings for a terradish graph
#'
#' Times the linear-system setup and solve used internally by
#' \code{\link{terradish_algorithm}} for a fixed graph and conductance vector.
#' This is intended as a lightweight hardware-specific diagnostic for choosing
#' direct factorization settings on large rasters.
#'
#' @param data A \code{terradish_graph} object.
#' @param conductance Optional numeric conductance vector with one value per
#'   graph vertex. Defaults to constant conductance.
#' @param factorization Direct solver factorization modes to compare. Supported
#'   values are \code{"simplicial_ldl"}, \code{"simplicial_ll"}, and
#'   \code{"supernodal_ll"}.
#' @param solve_backend Direct-solver RHS solve backend to compare.
#'   \code{"matrix"} uses Matrix's standard \code{solve()} method. The
#'   experimental \code{"cholmod_cpp"} backend calls Matrix's CHOLMOD C API
#'   directly from compiled code after the factorization has been built. The
#'   experimental \code{"cholmod_cpp_cached"} backend also keeps the CHOLMOD
#'   factorization/update state inside a compiled external pointer.
#' @param n_replicates Number of times to repeat each setting.
#' @param perm Should CHOLMOD use a fill-reducing permutation?
#'
#' @details
#' The benchmark isolates the sparse linear algebra step. It does not run a
#' full likelihood optimization, line search, measurement-model subproblem, or
#' Gaussian smoothing update. Use it to learn whether \code{"simplicial_ldl"},
#' \code{"simplicial_ll"}, or \code{"supernodal_ll"} is fastest on a particular
#' machine and graph size, then pass the preferred choice through
#' \code{solver_control = list(factorization = "...")} in \code{\link{terradish}}
#' when using the direct solver. Use \code{solve_backend} to evaluate whether
#' either experimental C++ CHOLMOD path is faster for the right-hand-side count
#' used by the graph. The cached backend is most informative when compared
#' inside a full optimization, where nearby evaluations can reuse the same
#' symbolic analysis.
#'
#' Very small rasters can give noisy timings because setup and R overhead are a
#' large fraction of elapsed time. The comparison is most informative on graphs
#' large enough that Cholesky factorization is a visible part of the fit time.
#'
#' @return A data frame with one row per benchmark replicate and columns for
#'   graph size, factorization mode, solve backend, setup time, solve time, and
#'   total time.
#'
#' @examples
#' library(terra)
#' data(melip)
#' melip.altitude <- terra::unwrap(melip.altitude)
#' melip.coords <- terra::unwrap(melip.coords)
#' surface <- conductance_surface(melip.altitude, melip.coords, directions = 8)
#' terradish_solver_benchmark(
#'   surface,
#'   factorization = c("simplicial_ldl", "simplicial_ll"),
#'   n_replicates = 1
#' )
#'
#' @export
terradish_solver_benchmark <- function(data,
                                       conductance = NULL,
                                       factorization = c("simplicial_ldl",
                                                         "simplicial_ll",
                                                         "supernodal_ll"),
                                       solve_backend = "matrix",
                                       n_replicates = 1L,
                                       perm = TRUE)
{
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  factorization <- match.arg(factorization, several.ok = TRUE)
  solve_backend <- match.arg(solve_backend,
                             c("matrix", "cholmod_cpp", "cholmod_cpp_cached"),
                             several.ok = TRUE)
  n_replicates <- as.integer(n_replicates)[1]
  if (is.na(n_replicates) || n_replicates < 1L)
    stop("`n_replicates` must be a positive integer", call. = FALSE)

  n_vertices <- nrow(data$x)
  if (is.null(conductance))
    conductance <- rep(1, n_vertices)
  conductance <- as.numeric(conductance)
  if (length(conductance) != n_vertices ||
      any(!is.finite(conductance)) ||
      any(conductance <= 0))
    stop("`conductance` must contain one finite positive value per graph vertex",
         call. = FALSE)

  rhs <- .graph_rhs(data, n_vertices)
  rows <- vector("list", length(factorization) * length(solve_backend) * n_replicates)
  row_id <- 0L

  for (mode in factorization)
  {
    for (backend in solve_backend)
    {
      for (replicate in seq_len(n_replicates))
      {
        elapsed_start <- proc.time()[["elapsed"]]
        solver_state <- .terradish_solver_setup(
          data,
          conductance,
          solver = "direct",
          solver_control = list(
            factorization = mode,
            solve_backend = backend,
            perm = isTRUE(perm)
          )
        )
        solve_result <- .terradish_solver_solve(solver_state, rhs)
        total_time <- proc.time()[["elapsed"]] - elapsed_start

        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          solver = "direct",
          factorization = mode,
          solve_backend = backend,
          replicate = replicate,
          n_vertices = n_vertices,
          n_rhs = ncol(rhs),
          setup_time = unname(solve_result$info$setup_time),
          solve_time = unname(solve_result$info$solve_time),
          total_time = unname(total_time),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows)
}
