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

# One exact, sparse-preserving Schur step: eliminate the local vertices `elim`
# from the (general, numerically symmetric) sparse operator `L`, returning the
# reduced operator on the kept vertices and the size of the interface that
# filled in. `L` is kept as a general dgCMatrix, not symmetric-packed, so the
# asymmetric tile subset `L[elim, keep]` stays sparse instead of forcing a dense
# coercion. Eliminating `elim` couples only the kept vertices adjacent to `elim`
# (the interface), so the dense work is `interface x interface`, not
# `keep x keep`. `L[elim, elim]` is a proper principal sub-Laplacian of a
# connected graph and so is positive definite (a Dirichlet block).
.schur_step_sparse <- function(L, elim)
{
  n    <- nrow(L)
  keep <- setdiff(seq_len(n), elim)
  Lee  <- forceSymmetric(L[elim, elim, drop = FALSE])
  Lek  <- as(L[elim, keep, drop = FALSE], "CsparseMatrix")
  Lkk  <- L[keep, keep, drop = FALSE]
  iface <- which(diff(Lek@p) > 0L)        # kept vertices coupled to `elim`
  m <- length(iface)
  if (m)
  {
    Lei <- Lek[, iface, drop = FALSE]
    Cee <- tryCatch(
      Cholesky(Lee, LDL = TRUE, perm = TRUE),
      error = function(e)
        stop("a tile interior block is not positive definite; tiles must ",
             "induce interior subgraphs connected to the retained set",
             call. = FALSE))
    U <- as.matrix(crossprod(Lei, solve(Cee, Lei)))   # interface x interface
    update <- sparseMatrix(i = rep.int(iface, m), j = rep(iface, each = m),
                           x = as.numeric(U), dims = dim(Lkk))
    Lkk <- Lkk - update
  }
  list(L = as(Lkk, "CsparseMatrix"), keep = keep, interface = m)
}

# Equal-population binning of a coordinate into `k` bins (quantile breaks).
.kron_bin <- function(z, k)
{
  if (k <= 1L) return(rep.int(1L, length(z)))
  br <- stats::quantile(z, probs = seq(0, 1, length.out = k + 1L), names = FALSE)
  br <- unique(br); br[1L] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(z, breaks = br, include.lowest = TRUE))
}

# Build a partition of all graph vertices into tiles. Honors a user-supplied
# `tiles` (a list of index vectors, or a per-vertex label vector); otherwise
# splits the vertices into a near-square grid of `n_tiles` spatial blocks from
# `vertex_coordinates`, falling back to contiguous index blocks.
.terradish_tile_groups <- function(data, tiles, n_tiles, n_vertices)
{
  if (!is.null(tiles))
  {
    if (is.list(tiles))
      groups <- lapply(tiles, as.integer)
    else
    {
      if (length(tiles) != n_vertices)
        stop("`tiles` label vector must have one entry per graph vertex",
             call. = FALSE)
      groups <- unname(split(seq_len(n_vertices), as.factor(tiles)))
    }
  }
  else
  {
    n_tiles <- max(1L, as.integer(n_tiles))
    coords  <- data$vertex_coordinates
    if (!is.null(coords) && nrow(coords) == n_vertices)
    {
      nb <- max(1L, floor(sqrt(n_tiles)))
      key <- interaction(.kron_bin(coords[, 1L], nb),
                         .kron_bin(coords[, 2L], nb), drop = TRUE)
      groups <- unname(split(seq_len(n_vertices), key))
    }
    else
    {
      width <- ceiling(n_vertices / n_tiles)
      blk   <- ceiling(seq_len(n_vertices) / width)
      groups <- unname(split(seq_len(n_vertices), blk))
    }
  }
  groups <- groups[vapply(groups, length, integer(1)) > 0L]
  if (!identical(sort(unlist(groups, use.names = FALSE)), seq_len(n_vertices)))
    stop("`tiles` must be a partition of all graph vertices", call. = FALSE)
  groups
}

#' Exact tiled (out-of-core) Schur/Kron reduction of a terradish graph
#'
#' Reduces a full terradish graph Laplacian onto a retained set of boundary
#' vertices (usually the focal sampling sites) by eliminating every other vertex
#' with an exact Schur complement, but does so tile by tile instead of in a
#' single interior factorization. Sequential Schur complements compose, so the
#' reduced Laplacian is identical to the one from
#' \code{\link{terradish_kron_reduce}}; the difference is purely computational.
#' Each tile contributes one small, local interior solve whose fill stays on
#' that tile's interface, so the peak working memory is one tile's interior
#' factor plus the (sparse) interface operator, never the full interior
#' factorization. That is the structure that lets the reduction proceed when the
#' full set of focal potentials no longer fits in memory: the regime past the
#' direct solver's memory wall and past the point where even an algebraic
#' multigrid solve must hold all focal potentials at once.
#'
#' @param data A \code{terradish_graph} object returned by
#'   \code{\link{conductance_surface}}.
#' @param conductance Numeric vertex conductance values, or a conductance-model
#'   result list containing a \code{conductance} element.
#' @param boundary Integer vertex indices to retain in the reduced graph.
#'   Defaults to the unique focal vertices in \code{data$demes}.
#' @param tiles Optional partition of the graph vertices into tiles, either a
#'   list of integer index vectors or a length-\code{nrow(data$x)} vector of
#'   per-vertex tile labels. When \code{NULL} (the default) the vertices are
#'   split into roughly \code{n_tiles} spatial blocks from
#'   \code{data$vertex_coordinates}, falling back to contiguous index blocks if
#'   coordinates are unavailable.
#' @param n_tiles Approximate number of tiles to build when \code{tiles} is
#'   \code{NULL}. The vertices are divided into a near-square grid of this many
#'   spatial blocks. More tiles lower the peak memory but raise the total work.
#' @param covariance If \code{TRUE}, also return the generalized inverse of the
#'   reduced Laplacian on the retained vertices, that is the focal
#'   effective-resistance covariance.
#'
#' @details The elimination is exact for any partition: \code{tiles} changes
#' only how much fill, and therefore how much peak memory, the reduction incurs.
#' Spatially compact tiles keep interfaces small. Because focal (boundary)
#' vertices are never eliminated, focal effective resistances are preserved
#' exactly, in contrast to approximate node lumping. The \code{peak} element of
#' the result reports the largest per-tile interior solve, the largest
#' interface, and the largest working-operator nonzero count, a direct read on
#' the memory the reduction actually needed.
#'
#' This promotes the tiled reduction as a standalone, exact primitive. It is not
#' yet wired into \code{\link{terradish_algorithm}} as a solver, because the
#' likelihood gradient needs the adjoint of conductance through the Schur
#' complement; that integration is a separate step.
#'
#' @return A list of class \code{"terradish_kron_reduction"} with the reduced
#'   \code{laplacian}, retained \code{boundary} vertices, eliminated
#'   \code{interior} vertices, the \code{tiles} used, a \code{peak} list of
#'   memory diagnostics (\code{interior}: largest per-tile interior solve;
#'   \code{interface}: largest interface; \code{nnz}: largest working-operator
#'   nonzero count), and optionally \code{covariance}.
#'
#' @seealso \code{\link{terradish_kron_reduce}} for the single-shot reduction
#'   this matches exactly.
#' @export
terradish_kron_reduce_tiled <- function(data,
                                        conductance,
                                        boundary = unique(data$demes),
                                        tiles = NULL,
                                        n_tiles = 16L,
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
    out <- list(laplacian = reduced, boundary = boundary, interior = interior,
                tiles = list(), n_boundary = length(boundary), n_interior = 0L,
                peak = list(interior = 0L, interface = 0L, nnz = length(Q@x)),
                covariance = if (isTRUE(covariance)) ginv(as.matrix(reduced)) else NULL)
    class(out) <- "terradish_kron_reduction"
    return(out)
  }

  groups <- .terradish_tile_groups(data, tiles, n_tiles, n_vertices)

  present  <- seq_len(n_vertices)        # global ids currently in `cur`
  # general (not symmetric-packed) sparse so asymmetric tile subsets stay sparse
  cur      <- as(Q, "generalMatrix")
  peak_int <- peak_iface <- 0L
  peak_nnz <- length(cur@x)
  for (g in groups)
  {
    elim_local <- match(intersect(g, interior), present)
    elim_local <- elim_local[!is.na(elim_local)]
    if (!length(elim_local)) next
    step       <- .schur_step_sparse(cur, elim_local)
    cur        <- step$L
    present    <- present[step$keep]
    peak_int   <- max(peak_int, length(elim_local))
    peak_iface <- max(peak_iface, step$interface)
    peak_nnz   <- max(peak_nnz, length(cur@x))
  }
  # eliminate any interior vertex a partial partition left behind
  leftover <- match(intersect(present, interior), present)
  leftover <- leftover[!is.na(leftover)]
  if (length(leftover))
  {
    step    <- .schur_step_sparse(cur, leftover)
    cur     <- step$L
    present <- present[step$keep]
  }

  pos     <- match(boundary, present)    # align ascending `present` to `boundary`
  reduced <- forceSymmetric(cur[pos, pos, drop = FALSE])

  out <- list(laplacian = reduced, boundary = boundary, interior = interior,
              tiles = groups, n_boundary = length(boundary),
              n_interior = length(interior),
              peak = list(interior = peak_int, interface = peak_iface, nnz = peak_nnz),
              covariance = if (isTRUE(covariance)) ginv(as.matrix(reduced)) else NULL)
  class(out) <- "terradish_kron_reduction"
  out
}
