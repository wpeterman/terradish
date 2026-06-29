test_that("conductance field basis and GMRF precision are well-formed", {
  set.seed(1)
  coords <- as.matrix(expand.grid(x = 1:5, y = 1:5))
  fld <- terradish:::.build_conductance_field(coords, G = 4, eps = 1e-3)
  expect_equal(nrow(fld$Z), nrow(coords))
  expect_equal(ncol(fld$Z), fld$m)
  # each cell assigned to exactly one occupied coarse cell
  expect_true(all(rowSums(fld$Z) == 1))
  # Q is symmetric PD (Laplacian + ridge)
  expect_equal(fld$Q, t(fld$Q), tolerance = 1e-12)
  expect_gt(min(eigen(fld$Q, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("penalized (theta, u) gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  dat <- melip_fixture(keep = 1:10)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  n <- length(surface$demes)
  X <- assemble_model_matrix(~ altitude + forestcover, surface$x)
  fld <- terradish:::.build_conductance_field(surface$vertex_coordinates,
                                              G = 3, eps = 1e-3)
  D <- cbind(X, fld$Z)
  cm <- terradish:::.design_loglinear_model(D)
  uidx <- (ncol(X) + 1L):(ncol(X) + fld$m)

  E <- as.matrix(terradish_algorithm(
    cm, leastsquares, surface, S = diag(n),
    theta = c(0.3, -0.3, rnorm(fld$m, 0, 0.2)),
    objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
  set.seed(7)
  S <- rWishart(1, df = 200, Sigma = (E + 0.3 * diag(n)) / 200)[, , 1]

  par <- c(0.2, -0.2, rnorm(fld$m, 0, 0.2))
  pen <- function(p, grad)
    terradish:::.hierarchical_penalized(
      p, tau2 = 0.5, cm = cm, measurement_model = wishart_covariance,
      data = surface, S = S, nu = 200, uidx = uidx, Q = fld$Q,
      nonnegative = TRUE, solver = "direct", solver_control = NULL,
      phi_state = NULL, want_grad = grad)
  an <- pen(par, TRUE)$gradient
  nd <- numDeriv::grad(function(p) pen(p, FALSE)$objective, par)
  expect_equal(an, nd, tolerance = 1e-4)
})

test_that("terradish_hierarchical reduces toward terradish as tau2 -> 0", {
  skip_on_cran()
  dat <- melip_fixture(keep = 1:16)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  n <- length(surface$demes)
  X <- as.matrix(surface$x[, c("altitude", "forestcover")])
  cm0 <- terradish:::.design_loglinear_model(cbind(X))
  E <- as.matrix(terradish_algorithm(
    cm0, leastsquares, surface, S = diag(n), theta = c(0.4, -0.5),
    objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
  set.seed(3)
  S <- rWishart(1, df = 800, Sigma = (E + 0.2 * diag(n)) / 800)[, , 1]

  fit_h <- terradish_hierarchical(S ~ altitude + forestcover, data = surface,
                                  measurement_model = wishart_covariance, nu = 800,
                                  field_resolution = 4, tau2 = 1e-5, verbose = FALSE)
  fit_t <- terradish(S ~ altitude + forestcover, data = surface,
                     conductance_model = loglinear_conductance,
                     measurement_model = wishart_covariance, nu = 800,
                     leverage = FALSE,
                     control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))

  expect_lt(max(abs(fit_h$u)), 0.05)            # field collapses
  expect_equal(unname(fit_h$theta), unname(coef(fit_t)), tolerance = 0.05)
})

test_that("fixed tau2 hierarchical fits expose marginal log-likelihood", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 4)

  fit_h <- suppressWarnings(
    terradish_hierarchical(dat$melip.Fst ~ altitude, data = surface,
                           measurement_model = leastsquares,
                           field_resolution = 2L, tau2 = 1,
                           maxit = 20L, verbose = FALSE)
  )

  expect_true("logml" %in% names(fit_h))
  expect_true(is.finite(fit_h$logml))
  expect_equal(fit_h$logml, fit_h$logML)
  expect_true(is.finite(fit_h$loglik))
})

test_that("terradish_hierarchical recovers an unmapped feature into the field", {
  skip_on_cran()
  dat <- melip_fixture(keep = 1:24)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  n  <- length(surface$demes)
  Vc <- surface$vertex_coordinates
  X  <- as.matrix(surface$x[, c("altitude", "forestcover")])

  xr <- range(Vc[, 1]); yr <- range(Vc[, 2])
  blob <- 1.6 * exp(-(((Vc[, 1] - (xr[1] + 0.35 * diff(xr)))^2 +
                       (Vc[, 2] - (yr[1] + 0.6 * diff(yr)))^2) /
                      (2 * (0.18 * sqrt(diff(xr)^2 + diff(yr)^2))^2)))
  blob <- blob - mean(blob)
  cm_true <- terradish:::.design_loglinear_model(cbind(X, blob))
  E <- as.matrix(terradish_algorithm(
    cm_true, leastsquares, surface, S = diag(n), theta = c(0.5, -0.6, 1.0),
    objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
  set.seed(11)
  S <- rWishart(1, df = 1500, Sigma = (E + 0.25 * diag(n)) / 1500)[, , 1]

  fit_h <- terradish_hierarchical(S ~ altitude + forestcover, data = surface,
                                  measurement_model = wishart_covariance, nu = 1500,
                                  field_resolution = 5, tau2 = 0.3, verbose = FALSE)
  u_cells <- as.vector(fit_h$field$Z %*% fit_h$u)
  expect_gt(cor(u_cells, blob), 0.6)            # field localizes the unmapped feature

  # conductance_field returns a raster with the field
  r <- conductance_field(fit_h, surface, type = "field")
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 1L)

  # S3 methods run
  expect_output(print(fit_h), "Hierarchical conductance surface")
  expect_length(coef(fit_h), 2L)
})
