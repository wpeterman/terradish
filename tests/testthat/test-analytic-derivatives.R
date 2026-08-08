# Finite-difference validation of the analytic derivatives returned by
# terradish_algorithm(). These are the derivatives the optimizer follows and the
# curvature that summary()/confint() turn into standard errors, so they are
# checked against numDeriv for every measurement model the package ships.
#
# The gradient is exact by the envelope theorem even when a nuisance parameter
# sits on a constraint, but the profile Hessian is not: it needs the implicit
# derivative dphi/dE, which is zero for any nuisance parameter pinned at an
# active box constraint. The `wishart_covariates` test below is the regression
# test for that case.

deriv_fixture <- function(keep = 1:8)
{
  dat <- melip_fixture(keep)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  list(surface = surface,
       f = loglinear_conductance(~ altitude + forestcover, surface$x),
       S_dist = ifelse(dat$melip.Fst < 0, 0, dat$melip.Fst),
       n_sites = length(keep))
}

# Compare the analytic gradient and Hessian at `theta` against Richardson
# extrapolated finite differences of the objective.
expect_derivatives_match <- function(f, g, surface, S, theta, nu = NULL,
                                     tol = 1e-5, ...)
{
  objective <- function(par)
    terradish_algorithm(f, g, surface, S, par, nu = nu,
                        gradient = FALSE, hessian = FALSE, partial = FALSE,
                        ...)$objective

  analytic <- terradish_algorithm(f, g, surface, S, theta, nu = nu,
                                  partial = FALSE, ...)

  numeric_gradient <- numDeriv::grad(objective, theta, method = "Richardson")
  numeric_hessian <- numDeriv::hessian(objective, theta, method = "Richardson")

  expect_equal(as.numeric(analytic$gradient), numeric_gradient,
               tolerance = tol)
  expect_equal(as.matrix(analytic$hessian), numeric_hessian,
               tolerance = tol, ignore_attr = TRUE)
  invisible(analytic)
}

test_that("leastsquares gradient and Hessian match finite differences", {
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  expect_derivatives_match(fx$f, leastsquares, fx$surface, fx$S_dist,
                           c(-0.3, 0.3))
})

test_that("generalized_wishart gradient and Hessian match finite differences", {
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  expect_derivatives_match(fx$f, generalized_wishart, fx$surface, fx$S_dist,
                           c(-0.3, 0.3), nu = 1000)
})

test_that("mlpe gradient and Hessian match finite differences", {
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  expect_derivatives_match(fx$f, mlpe, fx$surface, fx$S_dist, c(-0.3, 0.3))
})

test_that("wishart_covariance gradient and Hessian match finite differences", {
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  S_cov <- simulate_covariance_response(
    theta = c(0.15, -0.2), formula = ~ altitude + forestcover,
    data = fx$surface, conductance_model = loglinear_conductance,
    tau = 0.7, sigma = 0.05, nu = 40, seed = 1)$covariance
  expect_derivatives_match(fx$f, wishart_covariance, fx$surface, S_cov,
                           c(0.15, -0.2), nu = 40)
})

test_that("the profile Hessian is correct when a nuisance parameter is pinned at a bound", {
  # Regression test. `wishart_covariates` constrains each kernel coefficient to
  # be non-negative. When one is estimated at exactly 0 it cannot respond to a
  # perturbation of E, so it must be excluded from the profile-likelihood
  # curvature correction. Including it inflated the Hessian by ~1.6% here and
  # biased every standard error derived from it.
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  S_cov <- simulate_covariance_response(
    theta = c(0.15, -0.2), formula = ~ altitude + forestcover,
    data = fx$surface, conductance_model = loglinear_conductance,
    tau = 0.7, sigma = 0.05, nu = 40, seed = 1)$covariance

  set.seed(11)
  site_env <- data.frame(env = rnorm(fx$n_sites))
  g <- wishart_covariates(site_env, model = "wishart_covariance")
  theta <- c(0.15, -0.2)

  fit <- terradish_algorithm(fx$f, g, fx$surface, S_cov, theta, nu = 40,
                             partial = FALSE)
  # the fixture is chosen so the kernel coefficient really is on its bound;
  # `.free_parameter_index()` uses a sqrt(.Machine$double.eps) tolerance, so
  # anything below that counts as pinned
  expect_lt(abs(unname(fit$phi[2])), sqrt(.Machine$double.eps))
  # tau itself is interior, so this is a partially-active constraint set rather
  # than the fully degenerate `boundary = TRUE` case
  expect_false(isTRUE(fit$boundary))

  expect_derivatives_match(fx$f, g, fx$surface, S_cov, theta, nu = 40)

  # with the pinned coefficient removed, the curvature must agree with the
  # unconstrained base model, which is what the fit reduces to at lambda = 0
  base <- terradish_algorithm(fx$f, wishart_covariance, fx$surface, S_cov,
                              theta, nu = 40, partial = FALSE)
  expect_equal(as.matrix(fit$hessian), as.matrix(base$hessian),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("free nuisance parameters still contribute to the profile Hessian", {
  # Complement to the test above: `wishart_drift_covariates` estimates an
  # interior drift slope, so the correction must NOT be dropped there.
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  S_cov <- simulate_covariance_response(
    theta = c(0.15, -0.2), formula = ~ altitude + forestcover,
    data = fx$surface, conductance_model = loglinear_conductance,
    tau = 0.7, sigma = 0.05, nu = 40, seed = 1)$covariance

  set.seed(11)
  site_env <- data.frame(env = rnorm(fx$n_sites))
  g <- wishart_drift_covariates(site_env, model = "wishart_covariance")
  theta <- c(0.15, -0.2)

  fit <- terradish_algorithm(fx$f, g, fx$surface, S_cov, theta, nu = 40,
                             partial = FALSE)
  expect_true(all(is.finite(fit$phi)))
  expect_derivatives_match(fx$f, g, fx$surface, S_cov, theta, nu = 40)

  # the drift term genuinely changes the curvature relative to the base model
  base <- terradish_algorithm(fx$f, wishart_covariance, fx$surface, S_cov,
                              theta, nu = 40, partial = FALSE)
  expect_false(isTRUE(all.equal(as.matrix(fit$hessian),
                                as.matrix(base$hessian))))
})

test_that(".free_parameter_index flags only parameters on an active bound", {
  free <- terradish:::.free_parameter_index(
    phi = c(0.8, 0, -3, 2),
    lower = c(0, 0, -Inf, -Inf),
    upper = c(Inf, Inf, Inf, 2))
  expect_equal(free, c(TRUE, FALSE, TRUE, FALSE))

  # no bounds supplied means every parameter is free
  expect_true(all(terradish:::.free_parameter_index(c(1, 2, 3))))
})

test_that(".constrained_inverse_hessian zeroes pinned rows and columns", {
  H <- matrix(c(4, 1, 0,
                1, 3, 1,
                0, 1, 2), 3, 3)
  inv <- terradish:::.constrained_inverse_hessian(H, c(TRUE, FALSE, TRUE))
  expect_equal(inv[2, ], c(0, 0, 0))
  expect_equal(inv[, 2], c(0, 0, 0))
  expect_equal(inv[c(1, 3), c(1, 3)], solve(H[c(1, 3), c(1, 3)]))

  # all free reduces to the ordinary inverse
  expect_equal(terradish:::.constrained_inverse_hessian(H, rep(TRUE, 3)),
               solve(H))
})

test_that("linear and spline conductance derivatives match finite differences", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")

  data(melip, package = "terradish", envir = environment())
  covariates <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
  names(covariates) <- c("altitude", "forestcover")
  coords <- terra::unwrap(melip.coords)[1:8]
  S_dist <- ifelse(melip.Fst[1:8, 1:8] < 0, 0, melip.Fst[1:8, 1:8])

  # identity-link conductance requires strictly positive values, so use the
  # 0-1 rescaling rather than the mean-zero standardization
  surface_01 <- conductance_surface(scale_to_0_1(covariates), coords,
                                    directions = 8)
  f_linear <- linear_conductance(~ altitude + forestcover, surface_01$x)
  expect_derivatives_match(f_linear, leastsquares, surface_01, S_dist,
                           c(0.8, 1.2))

  surface_z <- conductance_surface(scale_covariates(covariates), coords,
                                   directions = 8)
  f_spline <- smooth_loglinear_conductance(~ s(altitude, df = 3) + forestcover,
                                           surface_z$x)
  n_par <- length(attr(f_spline, "default"))
  expect_derivatives_match(f_spline, leastsquares, surface_z, S_dist,
                           seq(0.05, 0.2, length.out = n_par))
})

test_that("pairwise-covariate measurement models match finite differences", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  fx <- deriv_fixture()
  set.seed(11)
  site_env <- data.frame(env = rnorm(fx$n_sites))

  expect_derivatives_match(fx$f, mlpe_covariates(site_env), fx$surface,
                           fx$S_dist, c(-0.3, 0.3))

  pairs <- t(utils::combn(fx$n_sites, 2))
  subset_pairs <- pairs[seq(1, nrow(pairs), by = 2), , drop = FALSE]
  expect_derivatives_match(fx$f,
                           pair_subset_measurement_model(leastsquares,
                                                         subset_pairs),
                           fx$surface, fx$S_dist, c(-0.3, 0.3))
})
