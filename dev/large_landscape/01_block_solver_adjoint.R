# =============================================================================
# Prototype A: matrix-free block iterative solve + streamed adjoint gradient
# -----------------------------------------------------------------------------
# Tier-1 engine from SCALING_DESIGN.md. Demonstrates, on terradish conventions
# (edge weight w_e = c_i + c_j, last-vertex grounding, E = A G), that:
#   * the focal covariance E and a scalar loss can be solved with a block,
#     preconditioner-reused iterative solver (no Cholesky factor materialised);
#   * the FULL gradient dl/dtheta needs exactly ONE forward + ONE adjoint block
#     solve, independent of the number of conductance parameters p, via a
#     streamed per-edge contraction (no dense (N-1) x (N-1) object).
#
# The math is verified to 3e-9 against finite differences in verify.py.
#
# INTEGRATION POINTS (terradish):
#   * build_reduced_laplacian()  <-> assemble_reduced_laplacian() (src/radish.cpp)
#   * block_pcg_jacobi()         <-> pcg_reduced_laplacian()      (src/radish.cpp)
#                                    or amg_reduced_laplacian_solve() for AMG
#   * adjoint_gradient()         replaces the gradient portion of
#                                .terradish_algorithm_derivative_chunk(); note it
#                                is O(1) solves vs the current O(p) Hessian loop.
#   * warm_start argument        <-> solver_warm_start in terradish_algorithm()
#
# Depends: Matrix. Base R otherwise.
# =============================================================================

suppressMessages(library(Matrix))

## ---- graph construction (terradish convention: w_e = c_i + c_j) -------------

# edges: integer matrix, 2 columns, 1-based vertex ids, i != j (undirected).
# conductance: length-N positive vector. Grounding drops the last vertex.
build_reduced_laplacian <- function(conductance, edges) {
  N <- length(conductance)
  w <- conductance[edges[, 1]] + conductance[edges[, 2]]
  # full symmetric Laplacian via triplets, then ground (drop vertex N)
  i <- c(edges[, 1], edges[, 2])
  j <- c(edges[, 2], edges[, 1])
  L <- sparseMatrix(i = i, j = j, x = c(-w, -w), dims = c(N, N))
  diag(L) <- -rowSums(L)
  forceSymmetric(L[-N, -N, drop = FALSE])
}

## ---- RHS Z and selection A (both constant in theta) -------------------------
# Z[, k] = -1/N everywhere, +1 at the row of focal site k (if not the grounded
# vertex). A maps a reduced solution to the focal covariance: E = A %*% G.
make_Z <- function(N, demes) {
  m <- length(demes)
  Z <- matrix(-1 / N, nrow = N - 1L, ncol = m)
  keep <- demes < N
  Z[cbind(demes[keep], which(keep))] <- Z[cbind(demes[keep], which(keep))] + 1
  Z
}
make_A <- function(N, demes) {
  m <- length(demes)
  A <- matrix(-1 / N, nrow = m, ncol = N - 1L)
  keep <- demes < N
  A[cbind(which(keep), demes[keep])] <- A[cbind(which(keep), demes[keep])] + 1
  A
}

## ---- block Jacobi-PCG (matrix-free in spirit; here L is a sparse matvec) -----
# Solves L X = B for all columns of B together, sharing the diagonal
# preconditioner. Swap in amg_reduced_laplacian_solve() for the production path.
# warm_start: optional (N-1) x ncol(B) initial guess (cross-Newton-step reuse).
block_pcg_jacobi <- function(L, B, tol = 1e-10, maxit = 5000, warm_start = NULL) {
  d <- diag(L)
  Minv <- function(R) R / d
  X <- if (is.null(warm_start)) matrix(0, nrow(B), ncol(B)) else warm_start
  R <- B - as.matrix(L %*% X)
  Zr <- Minv(R)
  P <- Zr
  rz_old <- colSums(R * Zr)
  bnorm <- pmax(sqrt(colSums(B^2)), 1)
  for (it in seq_len(maxit)) {
    AP <- as.matrix(L %*% P)
    alpha <- rz_old / colSums(P * AP)
    X <- X + sweep(P, 2, alpha, "*")
    R <- R - sweep(AP, 2, alpha, "*")
    if (all(sqrt(colSums(R^2)) <= tol * bnorm)) break
    Zr <- Minv(R)
    rz_new <- colSums(R * Zr)
    beta <- rz_new / rz_old
    P <- Zr + sweep(P, 2, beta, "*")
    rz_old <- rz_new
  }
  list(X = X, iter = it, resid = sqrt(colSums(R^2)) / bnorm)
}

## ---- forward pass: theta -> conductance -> E --------------------------------
# cond_fun(theta) -> length-N conductance; dcond_fun(theta) -> N x p (dc/dtheta).
forward_E <- function(theta, edges, demes, cond_fun, solver = block_pcg_jacobi,
                      warm_start = NULL) {
  c_vec <- cond_fun(theta)
  N <- length(c_vec)
  L <- build_reduced_laplacian(c_vec, edges)
  Z <- make_Z(N, demes)
  A <- make_A(N, demes)
  sol <- solver(L, Z, warm_start = warm_start)
  G <- sol$X
  list(L = L, G = G, E = A %*% G, A = A, c_vec = c_vec, solver_info = sol)
}

## ---- adjoint gradient: ONE extra solve, streamed per-edge contraction --------
# loss_grad_E(E) must return dl/dE (m x m). Returns dl/dtheta and the loss inputs.
# The (N-1)x(N-1) object dl/dL = -Lambda G^T is NEVER formed; we contract it
# edge-by-edge (this is the part that should be tiled by raster block at scale).
adjoint_gradient <- function(fwd, edges, dcond_fun, theta, dl_dE,
                             solver = block_pcg_jacobi) {
  N <- length(fwd$c_vec)
  G <- fwd$G                      # (N-1) x m
  A <- fwd$A
  dl_dG <- crossprod(A, dl_dE)    # A^T (dl/dE) : (N-1) x m
  Lam <- solver(fwd$L, dl_dG)$X   # adjoint solve, (N-1) x m (reuses factor if cached)

  # Streamed contraction. For reduced indices a,b: (dl/dL)_{ab} = -(Lam G^T)_{ab}.
  # We need only entries touched by edges. Pad to length N (grounded row/col = 0).
  red <- N - 1L
  pad_rows <- function(M) rbind(M, 0)          # row N -> zeros (grounded)
  Gp   <- pad_rows(G)                          # N x m
  Lamp <- pad_rows(Lam)                        # N x m
  ei <- edges[, 1]; ej <- edges[, 2]
  # diagonal pieces dii = -<Lam_i, G_i>; off-diagonal dij = -<Lam_i, G_j>
  dot <- function(Ua, Ub) rowSums(Ua * Ub)
  dii <- -dot(Lamp, Gp)                         # length N
  dl_dw <- dii[ei] + dii[ej] -
    ( -dot(Lamp[ei, , drop = FALSE], Gp[ej, , drop = FALSE]) ) -
    ( -dot(Lamp[ej, , drop = FALSE], Gp[ei, , drop = FALSE]) )
  # accumulate edge sensitivities to vertices (w_e = c_i + c_j => dw/dc = 1 each)
  dl_dc <- numeric(N)
  dl_dc <- dl_dc + tabulate(ei, N) * 0          # init
  dl_dc <- `+`(dl_dc, rowsum_fast(dl_dw, ei, N))
  dl_dc <- `+`(dl_dc, rowsum_fast(dl_dw, ej, N))
  dcdt <- dcond_fun(theta)                       # N x p
  grad <- as.numeric(crossprod(dcdt, dl_dc))     # p
  list(grad = grad, dl_dc = dl_dc, Lambda = Lam)
}

# fast scatter-add of values v into bins idx (1..N)
rowsum_fast <- function(v, idx, N) {
  out <- numeric(N)
  tapply_sum <- tapply(v, idx, sum)
  out[as.integer(names(tapply_sum))] <- tapply_sum
  out
}

## ============================ self-test ======================================
if (sys.nframe() == 0L) {
  set.seed(1)
  nr <- 6L; nc <- 7L; N <- nr * nc
  cid <- function(r, c) (r - 1L) * nc + c
  E_ <- list()
  for (r in 1:nr) for (cc in 1:nc) {
    if (cc < nc) E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r, cc + 1L))
    if (r < nr)  E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc))
    if (r < nr && cc < nc) E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc + 1L))
    if (r < nr && cc > 1)  E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc - 1L))
  }
  edges <- do.call(rbind, E_)
  demes <- c(1L, 6L, 18L, 31L, 42L); m <- length(demes)
  X <- matrix(rnorm(N * 2), N, 2)
  cond_fun  <- function(th) as.numeric(exp(X %*% th))
  dcond_fun <- function(th) { cc <- cond_fun(th); cc * X }
  theta <- c(0.3, -0.4)

  # a quadratic mock measurement model on E
  Btar <- crossprod(matrix(rnorm(m * m), m, m))
  loss_fun    <- function(E) 0.5 * sum((E - Btar)^2)
  loss_grad_E <- function(E) (E - Btar)

  fwd <- forward_E(theta, edges, demes, cond_fun)
  l   <- loss_fun(fwd$E)
  ag  <- adjoint_gradient(fwd, edges, dcond_fun, theta, loss_grad_E(fwd$E))$grad

  # finite-difference check
  fd <- numeric(2); h <- 1e-6
  for (k in 1:2) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    fd[k] <- (loss_fun(forward_E(tp, edges, demes, cond_fun)$E) -
              loss_fun(forward_E(tm, edges, demes, cond_fun)$E)) / (2 * h)
  }
  # direct-solve cross-check of the block PCG
  L <- build_reduced_laplacian(cond_fun(theta), edges)
  Gd <- as.matrix(solve(L, make_Z(N, demes)))
  Gp <- block_pcg_jacobi(L, make_Z(N, demes))$X

  cat(sprintf("loss                : %.6f\n", l))
  cat(sprintf("adjoint grad        : %s\n", paste(sprintf('% .6f', ag), collapse = " ")))
  cat(sprintf("finite-diff grad    : %s\n", paste(sprintf('% .6f', fd), collapse = " ")))
  cat(sprintf("block-PCG vs direct : %.2e\n", max(abs(Gp - Gd))))
  cat("PASS if both errors are tiny (<1e-5 grad, <1e-8 solve).\n")
}
