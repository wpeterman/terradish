test_that("Gaussian reporting and fixed-scale profiles respect fitted bounds", {
  set.seed(24)
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                    crs = "EPSG:3857")
  terra::values(r) <- rnorm(100); names(r) <- "x"
  graph <- conductance_surface(r, terra::xyFromCell(r, seq(3, 98, 5)),
                                directions = 4, saveStack = TRUE)
  factory <- gaussian_smoothed_loglinear_conductance(graph, sigma_lower = .4,
    sigma_upper = 1.5)
  sim <- simulate_covariance_response(c(x = .7, sigma.x = .8), ~x, graph,
    conductance_model = factory, tau = 1, sigma = .1, nu = 10000, seed = 14)
  S <- sim$covariance
  fit <- terradish(S ~ x, graph, conductance_model = factory,
    measurement_model = wishart_covariance, nu = 10000, optimizer = "newton")
  scale <- gaussian_scale_summary(fit)
  expect_equal(scale$sigma_lower, .4)
  expect_equal(scale$sigma_upper, 1.5)
  expect_true(all(is.na(summary(fit)$ztable["sigma.x", c("z value", "Pr(>|z|)")])))
  expect_gte(confint(fit)["sigma.x", 1], .4)
  expect_lte(confint(fit)["sigma.x", 2], 1.5)
  expect_equal(gaussian_scale_summary(slim_terradish(fit)), scale)
  bounded <- fit; bounded$mle$theta["sigma.x"] <- .401
  expect_true(summary(bounded)$sigma_table$near_bound)
  profile <- gaussian_scale_profile(fit, "x", n = 5)
  expect_lte(profile$interval[1], .8)
  expect_gte(profile$interval[2], .8)
  expect_equal(unname(confint(profile$fit)["sigma.x", ]), unname(profile$interval))
  expect_error(gaussian_scale_profile(slim_terradish(fit), "x"), "closures")
})
