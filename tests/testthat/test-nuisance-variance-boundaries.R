variance_fixture <- function() {
  E <- diag(c(.7, 1.3, 2, 3)) + matrix(.15, 4, 4)
  list(E = E, S = .8 * E + .2 * diag(4),
    control = NewtonRaphsonControl(ctol = 1e-8, ftol = 1e-8, maxit = 200, verbose = FALSE))
}

test_that("tiny log-variance starts cannot hide an improving positive variance", {
  f <- variance_fixture()
  start <- wishart_covariance(f$E, f$S)
  start$phi <- c(tau = .9, sigma = -50)
  fit <- radish_subproblem(wishart_covariance, f$E, f$S, 40, phi = start, control = f$control)
  expect_equal(fit$convergence, 0L)
  expect_equal(unname(fit$phi[, 1]), c(.8, log(.2)), tolerance = 1e-5)
  expect_false(fit$fit$residual_variance_zero_limit)
})

test_that("a verified zero-variance limit retains graph structure", {
  f <- variance_fixture(); S <- .8 * f$E
  start <- wishart_covariance(f$E, S)
  start$phi <- c(tau = .9, sigma = -50)
  fit <- radish_subproblem(wishart_covariance, f$E, S, 40, phi = start, control = f$control)
  expect_equal(fit$convergence, 0L)
  expect_true(fit$fit$residual_variance_zero_limit)
  expect_true(fit$boundary)
  expect_false(fit$no_structure_boundary)
  expect_equal(unname(fit$phi["tau", 1]), .8, tolerance = 1e-7)
  expect_lt(exp(fit$phi["sigma", 1]), 1e-300)
  expect_equal(fit$loglikelihood, wishart_covariance(f$E, S,
    phi = c(.8, -750), nu = 40)$objective, tolerance = 1e-8)
  expect_true(all(fit$jacobian_phi(diag(4))["sigma", ] == 0))
})

test_that("restarts respect fixed nuisance values and explicit variance bounds", {
  f <- variance_fixture(); start <- wishart_covariance(f$E, f$S)
  start$phi <- c(tau = .7, sigma = -50)
  start$lower[1] <- start$upper[1] <- .7
  fit <- radish_subproblem(wishart_covariance, f$E, f$S, 40, phi = start, control = f$control)
  expect_equal(unname(fit$phi["tau", 1]), .7)
  expect_gt(exp(fit$phi["sigma", 1]), .1)
  bounded <- wishart_covariance(f$E, .8 * f$E)
  bounded$lower[2] <- log(.01); bounded$phi[2] <- log(.01)
  fit <- radish_subproblem(wishart_covariance, f$E, .8 * f$E, 40,
    phi = bounded, control = f$control)
  expect_gte(fit$phi["sigma", 1], log(.01))
  expect_false(fit$fit$residual_variance_zero_limit)
  fixed <- start; fixed$lower <- fixed$upper <- fixed$phi
  fit <- radish_subproblem(wishart_covariance, f$E, f$S, 40, phi = fixed, control = f$control)
  expect_equal(c(fit$phi), unname(fixed$phi))
  expect_equal(fit$iters, 0L)
})

test_that("a warm-start exception preserves caller constraints", {
  f <- variance_fixture()
  g <- function(E, S, phi, ...) {
    if (!missing(phi) && phi[2] < -40) stop("deliberate warm-start failure")
    if (missing(phi)) wishart_covariance(E, S, ...) else wishart_covariance(E, S, phi, ...)
  }
  start <- g(f$E, f$S); start$phi <- c(tau = .7, sigma = -50)
  start$lower[1] <- start$upper[1] <- .7
  fit <- radish_subproblem(g, f$E, f$S, 40, phi = start, control = f$control)
  expect_equal(unname(fit$phi["tau", 1]), .7)
  expect_equal(fit$convergence, 0L)
})

test_that("distance and environmental Wishart fits receive the same safeguards", {
  f <- variance_fixture()
  D <- outer(diag(f$S), diag(f$S), "+") - 2 * f$S
  g <- generalized_wishart
  start <- g(f$E, D); start$phi <- c(tau = .9, sigma = -50)
  fit <- radish_subproblem(g, f$E, D, 40, phi = start, control = f$control)
  expect_equal(fit$convergence, 0L)
  expect_equal(unname(fit$phi[, 1]), c(.8, log(.2)), tolerance = 1e-5)
  for (kind in c("wishart_covariance", "generalized_wishart")) {
    g <- wishart_covariates(data.frame(environment = c(0, 1, 2, 4)), model = kind)
    S <- if (kind == "wishart_covariance") f$S else D
    start <- g(f$E, S); start$phi["sigma"] <- -50
    kernel <- grep("^lambda_", names(start$phi))
    start$lower[kernel] <- start$upper[kernel] <- 0
    start$phi[kernel] <- 0
    fit <- radish_subproblem(g, f$E, S, 40, phi = start, control = f$control)
    expect_equal(fit$convergence, 0L)
    expect_equal(unname(fit$phi[kernel, 1]), 0)
    expect_equal(exp(unname(fit$phi["sigma", 1])), .2, tolerance = 1e-5)
  }
})

test_that("zero-limit metadata reaches summaries without false Wald precision", {
  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5, crs = "EPSG:3857")
  terra::values(r) <- 0; names(r) <- "x"
  graph <- conductance_surface(r, terra::xyFromCell(r, c(3, 8, 17, 23)), directions = 4)
  f <- loglinear_conductance(~1, graph$x)
  E <- terradish_algorithm(f, wishart_covariance, graph, diag(4), theta = 0,
    phi = list(phi = c(tau = 1, sigma = log(.2)), lower = c(1, log(.2)), upper = c(1, log(.2))), nu = 40,
    gradient = FALSE, hessian = FALSE)$covariance
  S <- .8 * E
  start <- wishart_covariance(E, S); start$phi <- c(tau = .9, sigma = -50)
  # The deterministic profile start now reaches this endpoint directly.
  expect_no_warning(fit <- terradish(S ~ 1, graph, measurement_model = wishart_covariance, nu = 40,
    leverage = FALSE))
  expect_equal(fit$convergence$code, 0L)
  expect_equal(fit$fit$subproblem$convergence, 0L)
  expect_true(fit$fit$residual_variance_zero_limit)
  expect_false(.no_structure_boundary(fit$fit))
  x <- summary(fit)
  expect_true(all(is.na(x$phi_table["sigma", 2:4])))
  expect_true(x$phi_at_bound[["sigma"]])
  expect_match(x$residual_variance_note, "numerical zero limit")
})

test_that("partially fixed rho is not reset by the warm-start safeguard", {
  f <- variance_fixture()
  g <- function(E, S, phi, nu, gradient = TRUE, hessian = TRUE, partial = TRUE, nonnegative = TRUE) {
    if (missing(phi)) return(list(phi = c(tau = 1, rho = 0), lower = c(0, -Inf), upper = c(Inf, Inf)))
    dimnames_h <- list(c("tau", "rho"), c("tau", "rho"))
    list(objective = sum((phi - c(1, 9))^2), gradient = matrix(2 * (phi - c(1, 9)), 2, 1),
      hessian = matrix(c(2, 0, 0, 2), 2, dimnames = dimnames_h), boundary = FALSE,
      gradient_E = matrix(0, 4, 4), partial_E = matrix(0, 16, 2), partial_S = matrix(0, 2, 16),
      jacobian_E = function(x) x * 0, jacobian_S = function(x) x * 0)
  }
  start <- g(f$E, f$S); start$phi <- c(tau = 2, rho = 9)
  start$lower[2] <- start$upper[2] <- 9
  fit <- radish_subproblem(g, f$E, f$S, 40, phi = start, control = f$control)
  expect_equal(unname(fit$phi["rho", 1]), 9)
  expect_equal(fit$convergence, 0L)
})
