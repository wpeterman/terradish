# =============================================================================
# Prototype C: tiled Schur-complement (Kron) reduction to the focal sites
# -----------------------------------------------------------------------------
# The structural method B's results pointed to. Unlike adaptive node-lumping
# (prototype B), which removes internal path resistance and biased the focal
# covariance E by ~25%, Schur-complement elimination of non-focal cells is
# EXACT: the effective resistance among focal sites is invariant under Kron
# reduction. The win is computational, not statistical: interior cells are
# eliminated tile-by-tile with small, independent (parallelisable, out-of-core)
# local solves, leaving a small interface/focal system.
#
# This file provides:
#   schur_eliminate(L, elim)   exact Schur complement of L removing `elim`
#   schur_tiled(L, focal, tiles)  eliminate interior in tiles (sequential Schur)
#   eff_resistance(L, focal, ref) effective resistance among focal nodes
#
# `terradish::terradish_kron_reduce()` already exposes Kron reduction; this
# prototype demonstrates the *tiled* / domain-decomposition form and validates
# exactness against the full solve.
#
# Depends: Matrix. Validated in test_schur.R.
# =============================================================================
suppressMessages(library(Matrix))

# full (ungrounded) graph Laplacian from undirected edges with weights w_e
full_laplacian <- function(N, edges, w) {
  L <- sparseMatrix(i = c(edges[, 1], edges[, 2]),
                    j = c(edges[, 2], edges[, 1]),
                    x = c(-w, -w), dims = c(N, N))
  diag(L) <- -rowSums(L)
  forceSymmetric(L)
}

# Exact Schur complement of L onto keep = setdiff(1:N, elim):
#   S = L_kk - L_ke L_ee^{-1} L_ek
# For a connected graph and a proper interior `elim`, L_ee is SPD (Dirichlet
# Laplacian) so the elimination is exact and stable.
schur_eliminate <- function(L, elim) {
  N <- nrow(L); keep <- setdiff(seq_len(N), elim)
  if (length(elim) == 0L) return(list(S = L, keep = keep))
  Lee <- L[elim, elim, drop = FALSE]
  Lek <- L[elim, keep, drop = FALSE]
  Lkk <- L[keep, keep, drop = FALSE]
  S <- Lkk - crossprod(Lek, solve(Lee, Lek))   # L_ke = t(L_ek); exact
  list(S = as.matrix(forceSymmetric((S + t(S)) / 2)), keep = keep)
}

# Tiled elimination: remove all non-focal interior nodes, but in `tiles` groups,
# eliminating one group at a time from the partially reduced operator. Schur
# complements compose, so the result equals a single-shot elimination -- this is
# what makes interior elimination local/parallel/out-of-core.
schur_tiled <- function(L, focal, tiles) {
  N <- nrow(L)
  idx <- seq_len(N)                 # current global ids present in the operator
  cur <- L
  for (tile in tiles) {
    elim_global <- setdiff(tile, focal)             # never eliminate focal cells
    elim_local <- match(elim_global, idx)
    elim_local <- elim_local[!is.na(elim_local)]
    if (!length(elim_local)) next
    red <- schur_eliminate(cur, elim_local)
    cur <- Matrix(red$S, sparse = TRUE)
    idx <- idx[red$keep]
  }
  list(S = as.matrix(cur), keep = idx)              # idx are the surviving (focal) ids
}

# Effective resistance among focal nodes from a Laplacian L (ground at `ref`).
# R_ij = M_ii + M_jj - 2 M_ij with M = (L grounded at ref)^{-1}; R_i,ref = M_ii.
eff_resistance <- function(L, focal, ref = focal[1]) {
  N <- nrow(L)
  g <- which(seq_len(N) == ref)
  Lr <- as.matrix(L)[-g, -g, drop = FALSE]
  M <- solve(Lr)
  pos <- match(focal, setdiff(seq_len(N), ref))     # focal positions in grounded space (NA for ref)
  m <- length(focal); R <- matrix(0, m, m)
  d <- function(p) if (is.na(p)) 0 else M[p, p]
  for (a in 1:m) for (b in 1:m) {
    pa <- pos[a]; pb <- pos[b]
    maa <- d(pa); mbb <- d(pb)
    mab <- if (is.na(pa) || is.na(pb)) 0 else M[pa, pb]
    R[a, b] <- maa + mbb - 2 * mab
  }
  R
}
