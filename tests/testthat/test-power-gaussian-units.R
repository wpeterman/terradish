test_that("Gaussian power simulations use map-unit smoothing widths", {
  skip_on_cran()

  set.seed(11)
  raster <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 400,
                        ymin = 0, ymax = 400, crs = "local")
  noise <- matrix(rnorm(40 * 40), 40, 40)
  kernel <- outer(-6:6, -6:6,
                  function(row, col) exp(-(row^2 + col^2) / 8))
  terra::values(raster) <- as.vector(terra::focal(
    terra::rast(noise), w = kernel / sum(kernel), fun = sum, na.rm = TRUE)[])
  names(raster) <- "x"

  # Each cell is 10 map units wide; the requested sigma is four cells wide.
  xy <- expand.grid(x = seq(45, 355, by = 62),
                    y = seq(45, 355, by = 77.5))
  sites <- terra::vect(xy, geom = c("x", "y"), crs = "local")
  surface <- conductance_surface(raster, sites, directions = 8,
                                 saveStack = TRUE)
  gaussian <- gaussian_smoothed_loglinear_conductance(
    surface, scale_vars = "x")
  theta <- attr(gaussian(~ x, surface$x), "default")
  theta["x"] <- -0.8
  theta["sigma.x"] <- 40

  power <- covariance_response_power(
    theta = theta, formula = ~ x, data = surface, sample_sizes = 30,
    strategies = "spacefill", conductance_model = gaussian,
    fit_models = list(truth = list(formula = ~ x,
                                   conductance_model = gaussian,
                                   optimizer = "bfgs")),
    tau = 1, sigma = 0.05, nu = 5000, nsim = 1, seed = 3,
    control = NewtonRaphsonControl(maxit = 60L, verbose = FALSE)
  )

  sigma_estimate <- power$parameter_results$estimate[
    power$parameter_results$parameter == "sigma.x"]
  expect_true(isTRUE(power$results$fit_ok[[1]]))
  expect_equal(sigma_estimate, 40, tolerance = 0.25)
  expect_gt(power$results$conductance_cor[[1]], 0.95)
})
