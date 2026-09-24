test_that("spline rebuilds preserve knots, centers, and fitted coefficients", {
  set.seed(181)
  r <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 12)
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "x"
  coords <- terra::xyFromCell(r, c(1, 12, 25, 41, 64, 86, 99, 112, 133, 144))
  surface <- conductance_surface(r, coords, directions = 4, saveStack = TRUE)
  formula <- ~ s(x, df = 3)
  f <- smooth_loglinear_conductance(formula, surface$x)
  info <- attr(f, "smooth_loglinear_info")
  expect_length(info$smooth_specs[[1]]$centers, 3)
  stored <- attr(f, "plot_factory")
  rebuilt <- stored(formula, surface$x[10:60, , drop = FALSE])
  theta <- c(0.3, -0.2, 0.1)
  expect_equal(rebuilt(theta)$conductance, f(theta)$conductance[10:60], tolerance = 1e-12)
  old_factory <- get("smooth_loglinear_conductance", asNamespace("terradish"))
  old <- old_factory(formula, surface$x)
  log_shift <- log(f(theta)$conductance) - log(old(theta)$conductance)
  expect_lt(diff(range(log_shift)), 1e-12)
  fit_stub <- list(submodels = list(f_internal = f))
  expect_identical(.resolve_plot_conductance_model(fit_stub, smooth_loglinear_conductance), stored)
  expect_identical(.resolve_plot_conductance_model(fit_stub, NULL), stored)
  ci <- stored(formula, data.frame(x = min(surface$x$x)))(theta)$confint(theta, diag(3))
  expect_gt(ci[1, 2] - ci[1, 1], 0)
})

test_that("simulated Gaussian theta uses reported map units", {
  set.seed(182)
  r <- terra::rast(nrows = 8, ncols = 10, xmin = 0, xmax = 1000, ymin = 0, ymax = 800,
                   crs = "EPSG:32617")
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "x"
  coords <- terra::xyFromCell(r, c(1, 10, 28, 53, 71, 80))
  surface <- conductance_surface(r, coords, directions = 4, saveStack = TRUE)
  factory <- gaussian_smoothed_loglinear_conductance(surface, sigma_upper = 120)
  theta <- c(x = 0.3, sigma.x = 100)
  sim <- simulate_covariance_response(theta, ~x, surface, factory, nu = 30, sigma = 0.2, seed = 1)
  f <- factory(~x, surface$x)
  internal <- .conductance_model_to_internal(theta, f)
  expected <- terradish_algorithm(f, leastsquares, surface, S = diag(6), theta = internal,
                                   objective = FALSE, gradient = FALSE, hessian = FALSE,
                                   partial = FALSE)$covariance
  expect_equal(sim$E, as.matrix(expected), tolerance = 1e-10)
  expect_equal(c(sim$theta), unname(theta))
})
