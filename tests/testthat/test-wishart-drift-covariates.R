test_that("wishart_drift_covariates intercept-only equals wishart_covariance", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Sigma <- 0.8 * E + 0.2 * diag(3)
  set.seed(1)
  S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25

  g0 <- wishart_drift_covariates(model = "wishart_covariance")
  start <- g0(E, S, nu = 25)
  expect_equal(names(start$phi), c("tau", "sigma"))
  expect_equal(start$lower, c(0, -Inf))

  phi <- c(tau = 0.8, sigma = log(0.2))
  drift <- g0(E, S, phi = phi, nu = 25)
  base <- wishart_covariance(E, S, phi = phi, nu = 25)
  expect_equal(drift$objective, base$objective, tolerance = 1e-10)
  expect_equal(drift$gradient, base$gradient, tolerance = 1e-10)
  expect_equal(drift$hessian, base$hessian, tolerance = 1e-10)
  expect_equal(drift$gradient_E, base$gradient_E, tolerance = 1e-10)
  expect_equal(drift$partial_E, base$partial_E, tolerance = 1e-10)
  expect_equal(drift$partial_S, base$partial_S, tolerance = 1e-10)
})

test_that("wishart_drift_covariates with zero slope equals scalar nugget", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Sigma <- 0.8 * E + 0.2 * diag(3)
  set.seed(2)
  S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25

  g <- wishart_drift_covariates(data.frame(dens = c(-1, 0, 1)),
                                model = "wishart_covariance")
  base <- wishart_covariance(E, S, phi = c(tau = 0.8, sigma = log(0.2)), nu = 25)
  drift <- g(E, S, phi = c(tau = 0.8, sigma = log(0.2), gamma_dens = 0), nu = 25)
  expect_equal(drift$objective, base$objective, tolerance = 1e-10)
  expect_equal(drift$gradient_E, base$gradient_E, tolerance = 1e-10)
})

test_that("wishart_drift_covariates returns the covariance interface", {
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  g <- wishart_drift_covariates(data.frame(dens = c(-1, 0, 1)),
                                model = "wishart_covariance")
  Sigma <- 0.8 * E + diag(exp(c(-1.2, -1.0, -0.8)))
  set.seed(3)
  S <- rWishart(1, df = 30, Sigma = Sigma)[, , 1] / 30

  start <- g(E, S, nu = 30)
  fit <- g(E, S, phi = start$phi, nu = 30)

  expect_s3_class(g, "terradish_measurement_model")
  expect_equal(names(start$phi), c("tau", "sigma", "gamma_dens"))
  expect_equal(start$lower, c(0, -Inf, -Inf))
  expect_true(is.finite(fit$objective))
  expect_equal(dim(fit$gradient), c(3, 1))
  expect_equal(dim(fit$hessian), c(3, 3))
  expect_equal(dim(fit$gradient_E), dim(E))
  expect_equal(dim(fit$partial_E), c(length(E), 3))
  expect_equal(dim(fit$partial_S), c(3, sum(lower.tri(S, diag = TRUE))))
  expect_equal(dim(fit$jacobian_E(diag(nrow(E)))), dim(E))
  expect_equal(dim(fit$jacobian_S(diag(nrow(E)))), dim(S))
})

test_that("wishart_drift_covariates gradient/hessian match numDeriv (covariance)", {
  skip_if_not_installed("numDeriv")
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Z <- c(-1.3, 0.4, 0.9)
  Sigma <- 0.7 * E + diag(exp(cbind(1, scale(Z, scale = FALSE)) %*% c(log(0.25), 0.35)))
  set.seed(4)
  S <- rWishart(1, df = 25, Sigma = Sigma)[, , 1] / 25

  g <- wishart_drift_covariates(data.frame(dens = Z), model = "wishart_covariance")
  phi <- c(tau = 0.7, sigma = log(0.25), gamma_dens = 0.35)
  fit <- g(E, S, phi = phi, nu = 25)
  f_obj <- function(p) g(E, S, phi = p, nu = 25,
                         gradient = FALSE, hessian = FALSE, partial = FALSE)$objective

  expect_equal(c(fit$gradient), numDeriv::grad(f_obj, phi), tolerance = 1e-5)
  expect_equal(fit$hessian, numDeriv::hessian(f_obj, phi), tolerance = 1e-4,
               ignore_attr = TRUE)

  # gradient_E vs numerical
  f_objE <- function(Evec) {
    Em <- matrix(Evec, 3, 3); Em <- (Em + t(Em)) / 2
    g(Em, S, phi = phi, nu = 25,
      gradient = FALSE, hessian = FALSE, partial = FALSE)$objective
  }
  gE_num <- matrix(numDeriv::grad(f_objE, c(E)), 3, 3)
  gE_num <- (gE_num + t(gE_num)) / 2
  expect_equal(fit$gradient_E, gE_num, tolerance = 1e-4)
})

test_that("wishart_drift_covariates gradient/hessian match numDeriv (generalized)", {
  skip_if_not_installed("numDeriv")
  E <- matrix(c(1.2, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 1.1), nrow = 3, byrow = TRUE)
  Z <- c(-1.3, 0.4, 0.9)
  Sigma <- 0.7 * E + diag(exp(cbind(1, scale(Z, scale = FALSE)) %*% c(log(0.25), 0.3)))
  Sd <- diag(Sigma) %*% t(rep(1, 3)) + rep(1, 3) %*% t(diag(Sigma)) - 2 * Sigma

  g <- wishart_drift_covariates(data.frame(dens = Z), model = "generalized_wishart")
  phi <- c(tau = 0.7, sigma = log(0.25), gamma_dens = 0.3)
  fit <- g(E, Sd, phi = phi, nu = 25)
  f_obj <- function(p) g(E, Sd, phi = p, nu = 25,
                         gradient = FALSE, hessian = FALSE, partial = FALSE)$objective

  expect_equal(c(fit$gradient), numDeriv::grad(f_obj, phi), tolerance = 1e-5)
  expect_equal(fit$hessian, numDeriv::hessian(f_obj, phi), tolerance = 1e-4,
               ignore_attr = TRUE)
})

test_that("wishart_drift_covariates supports site subsetting", {
  site_env <- data.frame(dens = c(-1, 0, 1, 2))
  g <- wishart_drift_covariates(site_env, model = "wishart_covariance")
  sub_g <- attr(g, "subsetter")(c(1, 3, 4))
  expect_s3_class(sub_g, "terradish_measurement_model")
  expect_equal(nrow(attr(sub_g, "drift_covariates")), 3)
})

test_that("terradish recovers a drift slope and prefers it over a scalar nugget", {
  skip_on_cran()
  dat <- melip_fixture(keep = 1:24)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  n <- length(surface$demes)
  theta_true <- c(altitude = 0.2, forestcover = -0.25)

  # true covariance with a per-site drift surface driven by a site covariate
  alg0 <- terradish_algorithm(
    loglinear_conductance(~ altitude + forestcover, surface$x),
    leastsquares, surface, S = diag(n), theta = theta_true,
    objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)
  E_true <- as.matrix(alg0$covariance)
  Zc <- scale(seq_len(n))[, 1]
  nug <- as.vector(exp(cbind(1, scale(Zc, scale = FALSE)) %*% c(log(0.25), 0.7)))
  nu <- 500
  set.seed(11)
  S <- rWishart(1, df = nu, Sigma = (0.8 * E_true + diag(nug)) / nu)[, , 1]

  g <- wishart_drift_covariates(Zc, model = "wishart_covariance")
  fit_d <- suppressWarnings(terradish(
    S ~ altitude + forestcover, data = surface,
    conductance_model = loglinear_conductance, measurement_model = g,
    nu = nu, leverage = FALSE,
    control = NewtonRaphsonControl(maxit = 100, verbose = FALSE)))
  fit_s <- suppressWarnings(terradish(
    S ~ altitude + forestcover, data = surface,
    conductance_model = loglinear_conductance, measurement_model = wishart_covariance,
    nu = nu, leverage = FALSE,
    control = NewtonRaphsonControl(maxit = 100, verbose = FALSE)))

  expect_true("gamma_var1" %in% rownames(fit_d$fit$phi))
  expect_equal(unname(fit_d$fit$phi["gamma_var1", 1]), 0.7, tolerance = 0.2)
  aic <- function(f) -2 * f$loglik + 2 * (length(f$mle$theta) + nrow(f$fit$phi))
  expect_lt(aic(fit_d), aic(fit_s))
})
