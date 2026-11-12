# Every measurement model called without `phi` returns named starting values.
# generalized_wishart() used to return an unnamed vector, so fits with that
# model stored unnamed nuisance estimates.

test_that("all measurement models return named starting values", {
  # A small positive definite graph covariance and a noisy response, so the
  # least-squares starting regression is not an exact fit.
  set.seed(11)
  A <- matrix(rnorm(36), 6, 6)
  E <- crossprod(A) / 6 + diag(6)
  S_cov <- rWishart(1, df = 30, Sigma = 0.8 * E + 0.2 * diag(6))[, , 1] / 30
  S_dist <- dist_from_cov(S_cov)

  expect_named(generalized_wishart(E, S_dist, nu = 20)$phi, c("tau", "sigma"))
  expect_named(wishart_covariance(E, S_cov, nu = 20)$phi, c("tau", "sigma"))
  expect_named(leastsquares(E, S_dist)$phi, c("alpha", "beta", "tau"))
  expect_named(mlpe(E, S_dist)$phi, c("alpha", "beta", "tau", "rho"))

  # The values and bounds of the generalized Wishart start are unchanged.
  start <- generalized_wishart(E, S_dist, nu = 20)
  expect_identical(unname(start$phi), c(1, 0))
  expect_identical(start$lower, c(0, -Inf))
  expect_identical(start$upper, c(Inf, Inf))
})

test_that("generalized Wishart fits label their nuisance estimates", {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
  terra::values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
  names(r) <- "x"
  sites <- terra::xyFromCell(r, c(1, 8, 20, 37, 57, 64))
  graph <- conductance_surface(r, sites)
  simulated <- simulate_covariance_response(c(x = 0.6), ~x, graph,
                                            tau = 1, sigma = 0.1,
                                            nu = 1000, seed = 93)
  response <- dist_from_cov(simulated$covariance)

  fit <- terradish(response ~ x, graph,
                   measurement_model = generalized_wishart, nu = 1000)
  expect_identical(rownames(fit$fit$phi), c("tau", "sigma"))
  expect_identical(rownames(summary(fit)$phi_table), c("tau", "sigma"))

  # Naming the start does not change the likelihood: the same response as a
  # covariance gives the same estimates through wishart_covariance().
  covariance_fit <- terradish(simulated$covariance ~ x, graph,
                              measurement_model = wishart_covariance,
                              nu = 1000)
  expect_equal(unname(coef(fit)), unname(coef(covariance_fit)),
               tolerance = 1e-5)
  expect_equal(unname(fit$fit$phi[, 1]), unname(covariance_fit$fit$phi[, 1]),
               tolerance = 1e-5)
})
