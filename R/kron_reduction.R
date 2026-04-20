.terradish_full_laplacian <- function(data, conductance)
{
  n_vertices <- length(conductance)

  if (!is.null(data$laplacian))
  {
    Q <- data$laplacian
    Q@x[] <- -conductance[data$adj[1, ] + 1L] -
      conductance[data$adj[2, ] + 1L]
    Qd <- Diagonal(n_vertices, x = -rowSums(Q))
    return(forceSymmetric(Q + Qd))
  }

  if (is.null(data$adj))
    stop("`data` does not contain graph adjacency information", call. = FALSE)

  edge_i <- data$adj[1, ] + 1L
  edge_j <- data$adj[2, ] + 1L
  edge_conductance <- conductance[edge_i] + conductance[edge_j]
  Q <- sparseMatrix(
    i = c(edge_i, edge_j),
    j = c(edge_j, edge_i),
    x = -rep(edge_conductance, 2L),
    dims = c(n_vertices, n_vertices)
  )
  Q + Diagonal(n_vertices, x = -rowSums(Q))
}

#' Experimental Schur/Kron reduction of a terradish graph
#'
#' Reduces a full terradish graph Laplacian onto a retained set of boundary
#' vertices, usually the focal sampling sites, by eliminating all other graph
#' vertices with an exact Schur complement.
#'
#' @param data A \code{terradish_graph} object returned by
#'   \code{\link{conductance_surface}}.
#' @param conductance Numeric vertex conductance values, or a conductance-model
#'   result list containing a \code{conductance} element.
#' @param boundary Integer vertex indices to retain in the reduced graph.
#'   Defaults to the unique focal vertices in \code{data$demes}.
#' @param covariance If \code{TRUE}, also return the generalized inverse of the
#'   reduced Laplacian on the retained vertices.
#'
#' @details This is an experimental diagnostic/prototyping helper. The returned
#' Kron-reduced Laplacian is dense in general, so this is not automatically
#' faster than the main \code{\link{terradish_algorithm}} solve. Its purpose is
#' to make Schur/Kron reduction explicit so future approximations can be tested
#' against an exact reduction on small to moderate graphs.
#'
#' @return A list with the reduced \code{laplacian}, retained \code{boundary}
#'   vertices, eliminated \code{interior} vertices, and optionally
#'   \code{covariance}.
#' @export
terradish_kron_reduce <- function(data,
                                  conductance,
                                  boundary = unique(data$demes),
                                  covariance = FALSE)
{
  stopifnot(inherits(data, c("terradish_graph", "radish_graph")))
  if (is.list(conductance) && !is.null(conductance$conductance))
    conductance <- conductance$conductance
  conductance <- as.numeric(conductance)

  n_vertices <- nrow(data$x)
  if (length(conductance) != n_vertices)
    stop("`conductance` must have one value per graph vertex", call. = FALSE)

  boundary <- as.integer(boundary)
  if (length(boundary) < 2L || anyNA(boundary))
    stop("`boundary` must contain at least two retained vertex indices", call. = FALSE)
  if (any(boundary < 1L | boundary > n_vertices))
    stop("`boundary` values must be valid 1-based graph vertex indices", call. = FALSE)
  if (anyDuplicated(boundary))
    stop("`boundary` values must be unique", call. = FALSE)

  interior <- setdiff(seq_len(n_vertices), boundary)
  Q <- .terradish_full_laplacian(data, conductance)

  if (!length(interior))
  {
    reduced <- forceSymmetric(Q[boundary, boundary, drop = FALSE])
  }
  else
  {
    Qbb <- Q[boundary, boundary, drop = FALSE]
    Qbi <- Q[boundary, interior, drop = FALSE]
    Qib <- Q[interior, boundary, drop = FALSE]
    Qii <- forceSymmetric(Q[interior, interior, drop = FALSE])
    Qii_factor <- Cholesky(Qii, LDL = TRUE, perm = TRUE)
    reduced <- forceSymmetric(Qbb - Qbi %*% solve(Qii_factor, Qib))
  }

  out <- list(
    laplacian = reduced,
    boundary = boundary,
    interior = interior,
    n_boundary = length(boundary),
    n_interior = length(interior),
    covariance = if (isTRUE(covariance)) ginv(as.matrix(reduced)) else NULL
  )
  class(out) <- "terradish_kron_reduction"
  out
}
