# Small numerical checks run even when the larger solver/derivative suites skip.
cran_numeric_fixture <- function() {
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6,
                   ymin = 0, ymax = 6, crs = "EPSG:3857")
  terra::values(r) <- sin(seq_len(36) / 4) + seq_len(36) / 36
  names(r) <- "x"
  graph <- conductance_surface(r, terra::xyFromCell(r, c(1, 6, 13, 24, 31, 36)))
  S <- simulate_covariance_response(c(x = 0.4), ~x, graph,
    tau = 1, sigma = 0.1, nu = 100, seed = 96)$covariance
  list(graph = graph, S = S, f = loglinear_conductance(~x, graph$x))
}

test_that("small contrast-Wishart derivatives match central differences on CRAN", {
  a <- cran_numeric_fixture()
  evaluate <- function(theta, derivatives = TRUE) {
    terradish_algorithm(a$f, wishart_covariance, a$graph, a$S,
      theta = theta, nu = 100, gradient = derivatives,
      hessian = derivatives, partial = FALSE)
  }
  theta <- 0.2
  step <- 1e-4
  fit <- evaluate(theta)
  left <- evaluate(theta - step, FALSE)$objective
  right <- evaluate(theta + step, FALSE)$objective
  expect_equal(as.numeric(fit$gradient), (right - left) / (2 * step),
               tolerance = 1e-5)
  expect_equal(as.numeric(fit$hessian),
    (right - 2 * fit$objective + left) / step^2, tolerance = 2e-4)
})

test_that("small multilevel AMG solves agree with direct solves on CRAN", {
  a <- cran_numeric_fixture()
  direct <- terradish_algorithm(a$f, wishart_covariance, a$graph, a$S,
    theta = 0.2, nu = 100, partial = FALSE, solver = "direct")
  amg <- terradish_algorithm(a$f, wishart_covariance, a$graph, a$S,
    theta = 0.2, nu = 100, partial = FALSE, solver = "amg",
    solver_control = list(tol = 1e-10, maxit = 400L, coarse_enough = 8L))
  expect_identical(amg$solver_info$type, "amg")
  expect_true(all(amg$solver_info$converged))
  expect_equal(amg$objective, direct$objective, tolerance = 1e-7)
  expect_equal(as.matrix(amg$covariance), as.matrix(direct$covariance), tolerance = 1e-7)
  expect_equal(amg$gradient, direct$gradient, tolerance = 1e-6)
  expect_equal(amg$hessian, direct$hessian, tolerance = 1e-5)
})
