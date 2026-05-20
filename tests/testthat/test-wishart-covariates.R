test_that("wishart_covariates builds PSD kernels from site covariates", {
  site_env <- data.frame(altitude = c(-1, 0, 1, 2),
                         moisture = c(2, 1, 0, -1))
  g <- wishart_covariates(site_env, model = "wishart_covariance",
                          scale = TRUE)
  kernels <- attr(g, "kernel_covariates")

  expect_s3_class(g, "terradish_measurement_model")
  expect_equal(dim(kernels), c(4, 4, 2))
  expect_equal(dimnames(kernels)[[3]],
               c("kernel_altitude", "kernel_moisture"))
  for (k in seq_len(dim(kernels)[3]))
    expect_gte(min(eigen(kernels[, , k], symmetric = TRUE,
                         only.values = TRUE)$values), -1e-8)
})

test_that("wishart_covariates returns the covariance Wishart interface", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  site_env <- data.frame(altitude = c(-1, 0, 1))
  g <- wishart_covariates(site_env, model = "wishart_covariance")
  K <- attr(g, "kernel_covariates")[, , 1]
  Sigma <- 0.8 * E + 0.3 * K + 0.2 * diag(3)
  S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25

  start <- g(E, S, nu = 25)
  fit <- g(E, S, phi = start$phi, nu = 25)

  expect_equal(names(start$phi), c("tau", "lambda_altitude", "sigma"))
  expect_equal(start$lower, c(0, 0, -Inf))
  expect_true(is.finite(fit$objective))
  expect_equal(dim(fit$fitted), dim(S))
  expect_equal(dim(fit$gradient), c(3, 1))
  expect_equal(dim(fit$hessian), c(3, 3))
  expect_equal(dim(fit$gradient_E), dim(E))
  expect_equal(dim(fit$partial_E), c(length(E), 3))
  expect_equal(dim(fit$partial_S), c(3, sum(lower.tri(S, diag = TRUE))))
  expect_equal(dim(fit$jacobian_E(diag(nrow(E)))), dim(E))
  expect_equal(dim(fit$jacobian_S(diag(nrow(E)))), dim(S))
})

test_that("wishart_covariates matches base covariance Wishart with zero kernel weights", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Sigma <- 0.8 * E + 0.2 * diag(3)
  S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25
  g <- wishart_covariates(data.frame(altitude = c(-1, 0, 1)),
                          model = "wishart_covariance")

  wrapped <- g(E, S, phi = c(tau = 0.8, lambda_altitude = 0,
                             sigma = log(0.2)), nu = 25)
  base <- wishart_covariance(E, S, phi = c(tau = 0.8, sigma = log(0.2)),
                             nu = 25)

  expect_equal(wrapped$objective, base$objective, tolerance = 1e-8)
  expect_equal(wrapped$fitted, base$fitted, tolerance = 1e-8)
  expect_equal(wrapped$gradient_E, base$gradient_E, tolerance = 1e-8)
})

test_that("wishart_covariates supports generalized Wishart likelihoods", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Sigma <- 0.8 * E + 0.2 * diag(3)
  S <- diag(Sigma) %*% t(matrix(1, 3, 1)) +
    matrix(1, 3, 1) %*% t(diag(Sigma)) - 2 * Sigma
  g <- wishart_covariates(data.frame(altitude = c(-1, 0, 1)),
                          model = "generalized_wishart")

  wrapped <- g(E, S, phi = c(tau = 0.8, lambda_altitude = 0,
                             sigma = log(0.2)), nu = 25)
  base <- generalized_wishart(E, S, phi = c(tau = 0.8, sigma = log(0.2)),
                              nu = 25)

  expect_true(is.finite(wrapped$objective))
  expect_equal(wrapped$objective, base$objective, tolerance = 1e-8)
  expect_equal(wrapped$fitted, base$fitted, tolerance = 1e-8)
  expect_equal(dim(wrapped$gradient), c(3, 1))
  expect_equal(dim(wrapped$partial_S), c(3, sum(lower.tri(S))))
})

test_that("wishart_covariates supports site subsetting", {
  site_env <- data.frame(altitude = c(-1, 0, 1, 2),
                         moisture = c(2, 1, 0, -1))
  g <- wishart_covariates(site_env, model = "wishart_covariance")
  sub_g <- attr(g, "subsetter")(c(1, 3, 4))
  kernels <- attr(sub_g, "kernel_covariates")

  expect_s3_class(sub_g, "terradish_measurement_model")
  expect_equal(dim(kernels), c(3, 3, 2))
})

test_that("terradish can optimize covariance responses with Wishart covariates", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta_true <- c(altitude = 0.15, forestcover = -0.2)

  sim <- simulate_covariance_response(
    theta = theta_true,
    formula = ~ altitude + forestcover,
    data = surface,
    conductance_model = loglinear_conductance,
    tau = 0.7,
    sigma = 0.15,
    nu = 40,
    seed = 4
  )
  genetic_cov <- sim$covariance

  g <- wishart_covariates(dat$covariates[[1]], coords = dat$coords,
                          model = "wishart_covariance", scale = TRUE)
  fit <- suppressWarnings(
    terradish(genetic_cov ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = g,
              nu = 40,
              leverage = FALSE,
              control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  )

  expect_s3_class(fit, "terradish")
  expect_true(is.finite(fit$loglik))
  expect_true("lambda_altitude" %in% rownames(fit$fit$phi))
})

test_that("terradish can optimize distance responses with generalized Wishart covariates", {
  dat <- melip_fixture(keep = 1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  genetic_dist <- dat$melip.Fst

  g <- wishart_covariates(dat$covariates[[1]], coords = dat$coords,
                          model = "generalized_wishart", scale = TRUE)
  fit <- suppressWarnings(
    terradish(genetic_dist ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = g,
              nu = 1000,
              leverage = FALSE,
              control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  )

  expect_s3_class(fit, "terradish")
  expect_true(is.finite(fit$loglik))
  expect_false(is.complex(fit$loglik))
  expect_true("lambda_altitude" %in% rownames(fit$fit$phi))
})
