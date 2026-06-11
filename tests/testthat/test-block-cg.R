# Tests for the block conjugate-gradient solver (block_cg_reduced_laplacian and
# the solver = "block_cg" dispatch). Uses a small, well-conditioned grid graph,
# the regime block-CG (Jacobi-preconditioned) is designed for; large
# heterogeneous Laplacians are the AMG solver's regime.

small_grid <- function(nr = 25L, nc = 25L, seed = 1L) {
  cid <- function(r, c) (r - 1L) * nc + c
  e <- list()
  for (r in 1:nr) for (cc in 1:nc) {
    if (cc < nc) e[[length(e) + 1L]] <- c(cid(r, cc), cid(r, cc + 1L))
    if (r < nr)  e[[length(e) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc))
    if (r < nr && cc < nc) e[[length(e) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc + 1L))
    if (r < nr && cc > 1)  e[[length(e) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc - 1L))
  }
  ep <- do.call(rbind, e); storage.mode(ep) <- "integer"
  set.seed(seed); N <- nr * nc
  cond <- as.numeric(exp(matrix(rnorm(N * 2), N, 2) %*% c(0.4, -0.3)))
  demes <- as.integer(c(1, nc, N - nc + 1L, N, (N + 1L) %/% 2L, 7L))
  Z <- matrix(-1 / N, N - 1L, length(demes)); kd <- demes < N
  Z[cbind(demes[kd], which(kd))] <- Z[cbind(demes[kd], which(kd))] + 1
  list(ep = ep, cond = cond, demes = demes, Z = Z, N = N)
}

test_that("block_cg_reduced_laplacian matches the direct Cholesky solve", {
  g <- small_grid()
  L <- Matrix::forceSymmetric(assemble_reduced_laplacian(g$cond, g$ep))
  Xd <- as.matrix(Matrix::solve(Matrix::Cholesky(L), g$Z))
  bcg <- block_cg_reduced_laplacian(g$Z, g$cond, g$ep, tol = 1e-10, maxit = 5000)
  expect_true(all(bcg$converged))
  expect_lt(max(abs(bcg$solution - Xd)), 1e-6)
})

test_that("block_cg accepts a warm start and still converges to the same solution", {
  g <- small_grid()
  L <- Matrix::forceSymmetric(assemble_reduced_laplacian(g$cond, g$ep))
  Xd <- as.matrix(Matrix::solve(Matrix::Cholesky(L), g$Z))
  # warm start from a nearby parameter value's solution (the cross-step reuse
  # scenario): the initial residual block stays full rank, so block-CG is stable
  warm <- block_cg_reduced_laplacian(g$Z, g$cond * 1.15, g$ep, tol = 1e-10, maxit = 5000)$solution
  bcg <- block_cg_reduced_laplacian(g$Z, g$cond, g$ep, x0 = warm, tol = 1e-10, maxit = 5000)
  expect_true(all(bcg$converged))
  expect_lt(max(abs(bcg$solution - Xd)), 1e-6)
})

test_that("block_cg shares one Krylov space across right-hand sides", {
  # the block solve converges in far fewer iterations than the summed
  # single-RHS iteration count (the point of a block method)
  g <- small_grid()
  bcg <- block_cg_reduced_laplacian(g$Z, g$cond, g$ep, tol = 1e-10, maxit = 5000)
  single <- pcg_reduced_laplacian(g$Z, g$cond, g$ep, tol = 1e-10, maxit = 5000)
  expect_true(all(bcg$converged))
  expect_lt(bcg$iterations, sum(single$iterations))
})

test_that("solver = 'block_cg' dispatch matches the direct solve", {
  g <- small_grid()
  s <- list(edge_pairs = g$ep, demes = g$demes)
  st_b <- .terradish_solver_setup(s, g$cond, solver = "block_cg",
                                  solver_control = list(tol = 1e-10, maxit = 5000))
  expect_equal(st_b$type, "block_cg")
  sol_b <- .terradish_solver_solve(st_b, g$Z)$solution
  L <- Matrix::forceSymmetric(assemble_reduced_laplacian(g$cond, g$ep))
  Xd <- as.matrix(Matrix::solve(Matrix::Cholesky(L), g$Z))
  expect_lt(max(abs(as.matrix(sol_b) - Xd)), 1e-6)
})
