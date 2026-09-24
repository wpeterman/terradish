test_that("outer and joint Gaussian scales use the same map-unit kernel", {
  set.seed(311)
  r <- terra::rast(nrows = 12, ncols = 18, xmin = 0, xmax = 1800,
                   ymin = 0, ymax = 1200, crs = "EPSG:32617")
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "x"
  cells <- c(1, 18, 31, 50, 73, 90, 112, 135, 151, 173, 199, 216)
  coords <- terra::xyFromCell(r, cells)
  surface <- conductance_surface(r, coords, directions = 4, saveStack = TRUE)
  factory <- gaussian_smoothed_loglinear_conductance(surface,
                                                      sigma_lower = 50, sigma_upper = 180)
  sim <- simulate_covariance_response(c(x = 0.7, sigma.x = 100), ~x, surface,
                                       factory, tau = 1, sigma = 0.2, nu = 100, seed = 1)
  S <- sim$Sigma
  control <- NewtonRaphsonControl(maxit = 200)
  joint <- terradish(S ~x, surface, factory, wishart_covariance, nu = 100,
                     optimizer = "bfgs", control = control)
  outer <- terradish_scale_optim(S ~x, r, coords, scales = 100,
                                  lower = 50, upper = 180, directions = 4,
                                  measurement_model = wishart_covariance, nu = 100,
                                  control = control, tol = 0.01, maxit = 5)
  expect_equal(outer$scale_fun, "terradish")
  expect_lt(abs(outer$par[[1]] / coef(joint)[["sigma.x"]] - 1), 0.01)
  expect_equal(outer$fit$df, joint$df)
  expect_equal(outer$fit$aic, -2 * outer$fit$loglik + 2 * joint$df)
  legacy <- .resolve_scale_fun("multiScaleR")
  a <- do.call(legacy$fun, c(list(r, scale = 100), legacy$args))
  b <- multiScaleR::kernel_scale.raster(r, sigma = 100, kernel = "gaussian")
  expect_equal(terra::values(a), terra::values(b))
})
