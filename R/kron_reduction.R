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

# 1-based undirected edge list (m x 2) for the graph, used to find tile
# separators (cells with a neighbour in a different tile).
.kron_edges <- function(data)
{
  ep <- data$edge_pairs
  if (is.null(ep))
  {
    if (is.null(data$adj))
      stop("graph has no edge list (`edge_pairs`/`adj`)", call. = FALSE)
    ep <- t(data$adj) + 1L          # `adj` is 2 x m, 0-based
  }
  matrix(as.integer(ep), ncol = 2L)
}

# Local interior elimination for one tile: factor the tile's private interior
# block and return its dense Schur contribution onto its interface cells
# (positions within the global interface `B`). Independent across tiles, so this
# is the unit of parallel work; the payload is only the tile's own blocks.
.kron_tile_schur <- function(p)
{
  Cii <- tryCatch(
    Matrix::Cholesky(p$Lii, LDL = TRUE, perm = TRUE),
    error = function(e)
      stop("a tile interior block is not positive definite; tiles must induce ",
           "interior subgraphs connected to the interface", call. = FALSE))
  list(Bt = p$Bt,
       Ct = as.matrix(Matrix::crossprod(p$LiBt, Matrix::solve(Cii, p$LiBt))))
}

# Undirected edge list (m x 2, 1-based, i < j) from a sparse operator's nonzero
# pattern. Used at each nested-dissection level to bisect and to find the
# vertices straddling the cut (the separator).
.kron_operator_edges <- function(L)
{
  L <- as(L, "CsparseMatrix")
  j <- rep.int(seq_len(ncol(L)), diff(L@p)); i <- L@i + 1L
  k <- i < j
  cbind(i[k], j[k])
}

# Recursive nested dissection. Reduces operator `L` onto `keep` (local vertex
# indices) by bisecting the interior with a thin separator, recursively reducing
# each half onto its boundary, assembling the half contributions onto the
# interface (a multifrontal "extend-add"), and finally eliminating the
# separator. The two halves share no edges, so the result is the exact Schur
# complement on `keep`; bisecting keeps every single factorization bounded by
# the separator width rather than the whole interior or separator skeleton.
# Returns the reduced Laplacian on `keep` and the largest factorization
# dimension encountered (`maxfac`).
.kron_nd_reduce <- function(L, keep, coords, leaf, maxdepth = 40L, depth = 0L, par = 1L)
{
  n <- nrow(L); km <- logical(n); km[keep] <- TRUE; int <- which(!km)
  if (!length(int))
    return(list(laplacian = forceSymmetric(L[keep, keep, drop = FALSE]), maxfac = 0L))
  direct <- function() {
    s <- .schur_step_sparse(as(L, "generalMatrix"), int)
    P <- match(keep, (seq_len(n))[s$keep])
    list(laplacian = forceSymmetric(s$L[P, P, drop = FALSE]), maxfac = length(int))
  }
  if (length(int) <= leaf || depth >= maxdepth) return(direct())

  ax  <- if (diff(range(coords[int, 1])) >= diff(range(coords[int, 2]))) 1L else 2L
  med <- stats::median(coords[int, ax])
  side <- ifelse(coords[, ax] <= med, -1L, 1L)
  ed   <- .kron_operator_edges(L); cross <- side[ed[, 1]] != side[ed[, 2]]
  sepm <- logical(n); sepm[ed[cross, 1]] <- TRUE; sepm[ed[cross, 2]] <- TRUE
  Wm <- sepm & !km; Lh <- side < 0 & !km & !Wm; Rh <- side > 0 & !km & !Wm
  if (!any(Lh) || !any(Rh)) return(direct())          # no clean split -> direct

  M <- which(Wm | km); posM <- integer(n); posM[M] <- seq_along(M)
  M_op <- as(L[M, M, drop = FALSE], "generalMatrix"); maxfac <- 0L

  # The two halves share no edges, so they reduce independently. Recurse on each
  # (carrying its own boundary), returning the half's update onto that boundary;
  # the parent applies both extend-adds. The recursion forks on Unix when a
  # parallel budget remains (each child gets half the budget, so the number of
  # concurrent workers stays near `par`); the result is identical to sequential.
  half_fn <- function(hm) {
    hi <- which(hm)
    he <- ed[hm[ed[, 1]] | hm[ed[, 2]], , drop = FALSE]
    nb <- setdiff(unique(as.integer(he)), hi); nb <- nb[(Wm | km)[nb]]   # boundary in M
    hv <- c(hi, nb)
    rec <- .kron_nd_reduce(L[hv, hv, drop = FALSE], match(nb, hv),
                           coords[hv, , drop = FALSE], leaf, maxdepth, depth + 1L,
                           par %/% 2L)
    list(nb = nb, U = as.matrix(L[nb, nb, drop = FALSE]) - as.matrix(rec$laplacian),
         maxfac = rec$maxfac)
  }
  res <- if (par >= 2L && .Platform$OS.type == "unix")
    parallel::mclapply(list(Lh, Rh), half_fn, mc.cores = 2L, mc.preschedule = FALSE)
  else
    lapply(list(Lh, Rh), half_fn)
  if (any(vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))))
    stop("a parallel nested-dissection subtree failed; rerun with cores = 1 to see the error",
         call. = FALSE)
  for (r in res) {
    bp <- posM[r$nb]; M_op[bp, bp] <- M_op[bp, bp] - r$U
    maxfac <- max(maxfac, r$maxfac)
  }
  M_op <- as(forceSymmetric(M_op), "generalMatrix")
  Wpos <- posM[which(Wm)]
  if (length(Wpos)) {
    s <- .schur_step_sparse(M_op, Wpos)
    P <- match(keep, M[(seq_along(M))[s$keep]])
    list(laplacian = forceSymmetric(s$L[P, P, drop = FALSE]), maxfac = max(maxfac, length(Wpos)))
  } else {
    P <- match(keep, M)
    list(laplacian = forceSymmetric(M_op[P, P, drop = FALSE]), maxfac = maxfac)
  }
}

#' Exact nested-dissection (out-of-core, optionally parallel) Schur/Kron reduction of a terradish graph
#'
#' Reduces a full terradish graph Laplacian onto a retained set of boundary
#' vertices (usually the focal sampling sites) with an exact Schur complement,
#' bounding peak memory by never factorizing the whole interior at once. By
#' default it uses recursive nested dissection: the interior is bisected by a
#' thin separator into two independent halves, each reduced recursively, and the
#' separator is eliminated last, so the largest single factorization is bounded
#' by a separator width rather than by the interior or the separator skeleton.
#' The reduced Laplacian is identical to the one from
#' \code{\link{terradish_kron_reduce}}; the difference is purely computational.
#' This lets the reduction proceed at scales where holding the full set of focal
#' potentials no longer fits in memory: past the direct solver's memory wall and
#' past the point where even an algebraic multigrid solve must hold all focal
#' potentials at once.
#'
#' @param data A \code{terradish_graph} object returned by
#'   \code{\link{conductance_surface}}.
#' @param conductance Numeric vertex conductance values, or a conductance-model
#'   result list containing a \code{conductance} element.
#' @param boundary Integer vertex indices to retain in the reduced graph.
#'   Defaults to the unique focal vertices in \code{data$demes}.
#' @param tiles Partition control. \code{NULL} (the default) uses recursive
#'   nested dissection driven by \code{data$vertex_coordinates}. Supplying an
#'   explicit partition instead (a list of integer index vectors, or a
#'   length-\code{nrow(data$x)} vector of per-vertex tile labels) uses a flat
#'   two-level substructuring over those tiles, which is the parallel path (see
#'   \code{cores}).
#' @param n_tiles Leaf size for the default nested-dissection path: bisection
#'   stops when a subdomain interior has at most \code{n_tiles} vertices. Smaller
#'   leaves give a deeper recursion with smaller leaf factorizations; the largest
#'   factorization overall is set by the separator widths, not by this value.
#' @param covariance If \code{TRUE}, also return the generalized inverse of the
#'   reduced Laplacian on the retained vertices, that is the focal
#'   effective-resistance covariance.
#' @param cores Maximum number of parallel worker processes. \code{cores = 1}
#'   (the default) is fully sequential. With \code{cores > 1} the default
#'   nested-dissection path forks its independent recursive halves on Unix (the
#'   budget halves down the recursion, so roughly \code{cores} subtrees reduce at
#'   once); the explicit-partition path forks (Unix) or uses a socket cluster
#'   (Windows) across tiles. The result is identical to the sequential run.
#'   Forking shares memory, so the nested-dissection path falls back to
#'   sequential where forking is unavailable (for example on Windows).
#'
#' @details The default nested-dissection path bisects the interior by its
#' longer coordinate axis, marks the vertices straddling the cut as the
#' separator, recursively reduces each half onto its boundary, assembles the half
#' contributions onto the interface (a multifrontal extend-add), and eliminates
#' the separator last. Because the two halves share no edges the result is exact
#' for any bisection, and recursing keeps every factorization bounded by a
#' separator width, so the separator reduction no longer dominates memory at
#' large scale.
#'
#' The explicit-partition path (\code{tiles} supplied) eliminates each tile's
#' private interior independently onto a shared interface and then reduces that
#' interface in a single step; it parallelizes across tiles (see \code{cores}),
#' but its one interface factorization grows with the separator count, so it is
#' best at moderate scale. Either way focal vertices are never eliminated, so
#' focal effective resistances are preserved exactly, in contrast to approximate
#' node lumping.
#'
#' This promotes the reduction as a standalone, exact primitive. It is not yet
#' wired into \code{\link{terradish_algorithm}} as a solver, because the
#' likelihood gradient needs the adjoint of conductance through the Schur
#' complement; that integration is a separate step.
#'
#' @return A list of class \code{"terradish_kron_reduction"} with the reduced
#'   \code{laplacian}, retained \code{boundary} vertices, eliminated
#'   \code{interior} vertices, the \code{tiles} used (\code{NULL} for nested
#'   dissection), the \code{cores} used, and a \code{peak} list of memory
#'   diagnostics: \code{interior} (the largest single factorization), plus
#'   \code{separators} and \code{interface} on the explicit-partition path, and
#'   \code{nnz}. Optionally \code{covariance}.
#'
#' @seealso \code{\link{terradish_kron_reduce}} for the single-shot reduction
#'   this matches exactly.
#' @export
terradish_kron_reduce_tiled <- function(data,
                                        conductance,
                                        boundary = unique(data$demes),
                                        tiles = NULL,
                                        n_tiles = 2000L,
                                        covariance = FALSE,
                                        cores = 1L)
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

  # Default path: recursive nested dissection. Bisects the interior all the way
  # down, so every single factorization is bounded by a separator width rather
  # than by the whole interior or the whole separator skeleton. `n_tiles` is the
  # leaf size: bisection stops when a subdomain interior has <= n_tiles vertices.
  if (is.null(tiles))
  {
    coords <- data$vertex_coordinates
    if (is.null(coords) || nrow(coords) != n_vertices)
      stop("nested dissection needs `data$vertex_coordinates`; pass an explicit ",
           "`tiles` partition to use the flat substructuring path instead", call. = FALSE)
    leaf <- max(64L, as.integer(n_tiles))
    nd <- .kron_nd_reduce(as(Q, "generalMatrix"), boundary, as.matrix(coords), leaf,
                          par = max(1L, as.integer(cores)))
    reduced <- nd$laplacian
    out <- list(laplacian = reduced, boundary = boundary, interior = interior,
                tiles = NULL, n_boundary = length(boundary),
                n_interior = length(interior), cores = max(1L, as.integer(cores)),
                peak = list(interior = nd$maxfac, separators = NA_integer_,
                            interface = NA_integer_, nnz = length(Q@x)),
                covariance = if (isTRUE(covariance)) ginv(as.matrix(reduced)) else NULL)
    class(out) <- "terradish_kron_reduction"
    return(out)
  }

  # Explicit partition: flat two-level substructuring (parallel over tiles).
  groups <- .terradish_tile_groups(data, tiles, n_tiles, n_vertices)
  cores  <- max(1L, as.integer(cores))

  # general (not symmetric-packed) sparse so asymmetric tile subsets stay sparse
  Qg <- as(Q, "generalMatrix")

  # tile id per vertex, then separators: a vertex with a neighbour in another
  # tile. The interface B is the separators plus the retained (focal) vertices.
  tile_of <- integer(n_vertices)
  for (t in seq_along(groups)) tile_of[groups[[t]]] <- t
  ep    <- .kron_edges(data)
  cross <- tile_of[ep[, 1]] != tile_of[ep[, 2]]
  Bmask <- logical(n_vertices)
  Bmask[ep[cross, 1]] <- TRUE; Bmask[ep[cross, 2]] <- TRUE
  Bmask[boundary] <- TRUE
  B <- which(Bmask)

  # Per tile, extract only its private-interior blocks (bounded payload). The
  # private interiors of different tiles share no edges, so each tile's interior
  # eliminates independently -- this is the parallel unit.
  payloads <- lapply(groups, function(cells)
  {
    It <- cells[!Bmask[cells]]
    if (!length(It)) return(NULL)
    LiB <- as(Qg[It, B, drop = FALSE], "CsparseMatrix")   # interior x interface
    Bt  <- which(diff(LiB@p) > 0L)                        # interface cols touching It
    if (!length(Bt)) return(NULL)
    list(Lii = forceSymmetric(Qg[It, It, drop = FALSE]),
         LiBt = LiB[, Bt, drop = FALSE], Bt = Bt, ni = length(It))
  })
  payloads <- Filter(Negate(is.null), payloads)

  nw <- min(as.integer(cores), length(payloads))
  contribs <- if (nw > 1L && .Platform$OS.type == "unix")
  {
    # fork: workers share the parent's memory, so no per-tile serialization
    parallel::mclapply(payloads, .kron_tile_schur, mc.cores = nw)
  }
  else if (nw > 1L)
  {
    # Windows: socket cluster (each worker is sent only its own tile blocks).
    # Note that serializing those blocks has real overhead, so the parallel win
    # only materializes when the per-tile factorizations are large relative to
    # it; otherwise `cores = 1` is faster.
    cl <- parallel::makeCluster(nw)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, requireNamespace("Matrix", quietly = TRUE))
    parallel::parLapply(cl, payloads, .kron_tile_schur)
  }
  else
  {
    lapply(payloads, .kron_tile_schur)
  }
  if (any(vapply(contribs, function(x) inherits(x, "try-error") || is.null(x), logical(1))))
    stop("a parallel tile reduction failed; rerun with cores = 1 to see the error",
         call. = FALSE)

  # Assemble the interface operator S_B = Q[B,B] - sum_t (tile Schur contribution).
  S_B <- as(Qg[B, B, drop = FALSE], "generalMatrix")
  if (length(contribs))
  {
    # assemble all tile contributions in one pass (avoid growing vectors with c())
    ii <- unlist(lapply(contribs, function(cc) rep.int(cc$Bt, length(cc$Bt))), use.names = FALSE)
    jj <- unlist(lapply(contribs, function(cc) rep(cc$Bt, each = length(cc$Bt))), use.names = FALSE)
    xx <- unlist(lapply(contribs, function(cc) as.numeric(cc$Ct)), use.names = FALSE)
    S_B <- S_B - sparseMatrix(i = ii, j = jj, x = xx, dims = dim(S_B))
  }
  S_B <- as(forceSymmetric(S_B), "generalMatrix")

  # Reduce the interface onto the focal vertices: eliminate the separators.
  foc_in_B <- match(boundary, B)
  elim_B   <- setdiff(seq_along(B), foc_in_B)
  if (length(elim_B))
  {
    step    <- .schur_step_sparse(S_B, elim_B)
    kept_B  <- (seq_along(B))[step$keep]
    reduced <- step$L
  }
  else
  {
    kept_B <- seq_along(B); reduced <- S_B
  }
  pos     <- match(boundary, B[kept_B])  # align to `boundary` order
  reduced <- forceSymmetric(reduced[pos, pos, drop = FALSE])

  peak_int <- if (length(payloads))
    max(vapply(payloads, function(p) p$ni, integer(1))) else 0L

  out <- list(laplacian = reduced, boundary = boundary, interior = interior,
              tiles = groups, n_boundary = length(boundary),
              n_interior = length(interior), cores = cores,
              peak = list(interior = peak_int, separators = length(B) - length(boundary),
                          interface = length(B), nnz = max(length(Qg@x), length(S_B@x))),
              covariance = if (isTRUE(covariance)) ginv(as.matrix(reduced)) else NULL)
  class(out) <- "terradish_kron_reduction"
  out
}
