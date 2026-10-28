mlpe_mode_fixture <- function() {
  b <- readRDS(test_path("fixtures", "wishart-nuisance-mode.rds"))
  list(E = b$E, S = dist_from_cov(b$S))
}

test_that("MLPE profiling and derivatives are invariant to resistance units", {
  b <- mlpe_mode_fixture()
  fit <- function(E, phi = NULL) radish_subproblem(mlpe, E, b$S,
    nu = NULL, phi = phi,
    control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10))
  base <- fit(b$E)
  expect_equal(base$convergence, 0)
  expect_equal(as.numeric(base$loglikelihood), -1687.120991, tolerance = 1e-7)
  D <- diag(seq(.001, .002, length.out = nrow(b$E)))
  for (units in c(1e-12, 1e12)) {
    other <- fit(b$E * units)
    expect_equal(other$convergence, 0)
    expect_equal(other$loglikelihood, base$loglikelihood, tolerance = 1e-12)
    expect_equal(other$phi[2] * units, base$phi[2], tolerance = 1e-8)
    expect_equal(other$gradient * units, base$gradient, tolerance = 1e-7)
    expect_equal(other$jacobian_E(D * units) * units,
                 base$jacobian_E(D), tolerance = 1e-6)
    expect_true(.free_parameter_index(other$phi,
      c(-Inf, 0, -Inf, -Inf), rep(Inf, 4))[2])
  }
  warm <- fit(b$E, c(100, 1e-20, 10, -80))
  expect_equal(warm$phi, base$phi, tolerance = 1e-12)
  expect_equal(warm$loglikelihood, base$loglikelihood, tolerance = 1e-12)
})

test_that("nuisance inversion and active bounds respect parameter units", {
  H <- diag(c(1e-24, 1, 1e24))
  expect_equal(.safe_invert(H), diag(c(1e24, 1, 1e-24)), tolerance = 1e-12)
  expect_equal(.free_parameter_index(c(0, 1e-12, 1e12),
    lower = c(0, 0, 0), upper = rep(Inf, 3)), c(FALSE, TRUE, TRUE))
  b <- mlpe_mode_fixture()
  fixed <- c(alpha = 1, beta = 1e-12, tau = 0, rho = -2)
  fit <- radish_subproblem(mlpe, b$E, b$S, nu = NULL,
    phi = list(phi = fixed, lower = fixed, upper = fixed))
  expect_equal(c(fit$phi), unname(fixed))
  expect_equal(fit$loglikelihood, mlpe(b$E, b$S, fixed)$objective)
})

test_that("a zero Newton Hessian does not cause a singular inversion", {
  f <- function(x, gradient, hessian) list(objective = x[1]^4 / 4 - x[1],
    gradient = matrix(x[1]^3 - 1), hessian = matrix(3 * x[1]^2))
  fit <- BoxConstrainedNewton(0, f)
  expect_equal(fit$convergence, 0)
  expect_equal(as.numeric(fit$par), 1, tolerance = 1e-6)
})

test_that("reference changes preserve known path resistances on a stiff graph", {
  n <- 12L
  pairs <- cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
  graph <- list(x = data.frame(x = seq_len(n)), demes = seq_len(n),
    edge_pairs = pairs, adj = t(pairs - 1L), reduced_index = seq_len(n - 1L))
  class(graph) <- c("terradish_graph", "radish_graph")
  C <- rep(1, n); C[3:6] <- 1e10
  positions <- c(0, cumsum(1 / (head(C, -1L) + tail(C, -1L))))
  expected <- abs(outer(positions, positions, "-"))
  rhs <- .graph_rhs(graph, n)
  for (backend in c("matrix", "cholmod_cpp", "cholmod_cpp_cached")) {
    state <- .terradish_solver_setup(graph, C, "direct",
      list(solve_backend = backend))
    solution <- .terradish_solver_solve(state, rhs)$solution
    E <- graph_rhs_crossprod(graph$demes, n, as.matrix(solution))
    actual <- outer(diag(E), diag(E), "+") - E - t(E)
    expect_equal(actual, expected, tolerance = 1e-10)
    expect_equal(state$reference_vertex, 3L)
  }
  legacy <- .terradish_solver_setup(graph, C, "direct", list(reference = "last"))
  expect_error(.terradish_solver_solve(legacy, rhs), "residual tolerance")
  unchecked <- .terradish_solver_setup(graph, C, "direct",
    list(reference = "last", residual_tolerance = Inf))
  expect_gt(.terradish_solver_solve(unchecked, rhs)$info$relative_residual, 1e-8)
})
