# Tests for the `curvature` argument (exact vs Gauss-Newton/Fisher) added to
# terradish_algorithm() and terradish().

curvature_surface <- function(keep = 1:8) {
  dat <- melip_fixture(keep)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  S <- ifelse(dat$melip.Fst < 0, 0, dat$melip.Fst)
  f <- loglinear_conductance(~ altitude + forestcover, surface$x)
  list(surface = surface, S = S, f = f)
}

test_that("exact and gauss_newton curvature share objective and gradient", {
  cs <- curvature_surface()
  theta <- c(-0.3, 0.3)
  ex <- terradish_algorithm(cs$f, leastsquares, cs$surface, cs$S, theta,
                            partial = FALSE, curvature = "exact")
  gn <- terradish_algorithm(cs$f, leastsquares, cs$surface, cs$S, theta,
                            partial = FALSE, curvature = "gauss_newton")
  expect_equal(ex$objective, gn$objective)
  expect_equal(ex$gradient, gn$gradient)
  # exact path must be untouched by the new option (it only guards the two
  # residual-weighted second-derivative terms)
  expect_false(isTRUE(all.equal(ex$hessian, gn$hessian)))
})

test_that("gauss_newton curvature is symmetric everywhere", {
  cs <- curvature_surface()
  for (g in list(leastsquares, generalized_wishart)) {
    gn <- terradish_algorithm(cs$f, g, cs$surface, cs$S, c(-0.3, 0.3),
                              nu = 1000, partial = FALSE,
                              curvature = "gauss_newton")
    expect_equal(gn$hessian, t(gn$hessian), tolerance = 1e-10)
  }
})

test_that("gauss_newton curvature is positive semidefinite at the optimum", {
  # equals the Fisher information at the optimum, so PSD there (away from the
  # optimum the measurement model's observed E-curvature can be indefinite).
  for (g in list(leastsquares, generalized_wishart)) {
    res <- fit_fixture(keep = 1:8, measurement_model = g, curvature = "exact", nu = 1000)
    cs <- curvature_surface()
    gn <- terradish_algorithm(cs$f, g, cs$surface, cs$S, coef(res$fit),
                              nu = 1000, partial = FALSE,
                              curvature = "gauss_newton")
    ev <- eigen(gn$hessian, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev >= -1e-6 * max(abs(ev))))
  }
})

test_that("gauss_newton converges to the exact Hessian at a well-fitting optimum", {
  cs <- curvature_surface()
  # generate a response from the model itself (+ tiny noise) so residuals ~ 0;
  # there the dropped second-derivative terms vanish and GN == exact Hessian.
  theta <- c(-0.4, 0.4)
  E <- terradish_algorithm(cs$f, leastsquares, cs$surface, cs$S, theta,
                           objective = FALSE, gradient = FALSE, hessian = FALSE,
                           partial = FALSE)$covariance
  E <- as.matrix(E)
  R <- outer(diag(E), diag(E), "+") - 2 * E
  set.seed(1)
  Smodel <- 0.5 + 0.8 * R + 1e-4 * (function(M){M[] <- rnorm(length(M)); (M + t(M)) / 2})(R)
  diag(Smodel) <- 0
  ex <- terradish_algorithm(cs$f, leastsquares, cs$surface, Smodel, theta,
                            partial = FALSE, curvature = "exact")$hessian
  gn <- terradish_algorithm(cs$f, leastsquares, cs$surface, Smodel, theta,
                            partial = FALSE, curvature = "gauss_newton")$hessian
  expect_lt(norm(ex - gn, "F") / norm(ex, "F"), 0.05)
})

# wishart_covariance takes a covariance response (not an Fst distance), so it
# needs its own fixture: a response simulated from the model itself. Raising the
# Wishart degrees of freedom nu concentrates the draw toward its mean, shrinking
# the residual so the dropped second-derivative terms vanish (see the
# convergence test below).
curvature_cov_fixture <- function(keep = 1:8, theta = c(0.15, -0.2), nu = 40) {
  dat <- melip_fixture(keep)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  f <- loglinear_conductance(~ altitude + forestcover, surface$x)
  S <- simulate_covariance_response(theta = theta, formula = ~ altitude + forestcover,
                                    data = surface, conductance_model = loglinear_conductance,
                                    tau = 0.7, sigma = 0.05, nu = nu, seed = 1)$covariance
  list(surface = surface, S = S, f = f, theta = theta, nu = nu)
}

test_that("wishart_covariance: gauss_newton shares objective/gradient and is symmetric + PSD", {
  cf <- curvature_cov_fixture(nu = 40)
  ex <- terradish_algorithm(cf$f, wishart_covariance, cf$surface, cf$S, cf$theta,
                            nu = cf$nu, partial = FALSE, curvature = "exact")
  gn <- terradish_algorithm(cf$f, wishart_covariance, cf$surface, cf$S, cf$theta,
                            nu = cf$nu, partial = FALSE, curvature = "gauss_newton")
  # the toggle (in the shared derivative chunk) changes only the Hessian
  expect_equal(ex$objective, gn$objective)
  expect_equal(ex$gradient, gn$gradient)
  expect_false(isTRUE(all.equal(ex$hessian, gn$hessian)))
  # GN Hessian is a Fisher-information form: symmetric and PSD
  expect_equal(gn$hessian, t(gn$hessian), tolerance = 1e-10)
  ev <- eigen(gn$hessian, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-6 * max(abs(ev))))
})

test_that("wishart_covariance: gauss_newton converges to the exact Hessian as residuals vanish", {
  # concentrate the Wishart draw (large nu) so the response is ~ the model mean
  # at theta; there the residual-weighted terms drop out and GN == exact.
  cf <- curvature_cov_fixture(nu = 4000)
  ex <- terradish_algorithm(cf$f, wishart_covariance, cf$surface, cf$S, cf$theta,
                            nu = cf$nu, partial = FALSE, curvature = "exact")$hessian
  gn <- terradish_algorithm(cf$f, wishart_covariance, cf$surface, cf$S, cf$theta,
                            nu = cf$nu, partial = FALSE, curvature = "gauss_newton")$hessian
  expect_lt(norm(ex - gn, "F") / norm(ex, "F"), 0.05)
})

test_that("terradish() accepts curvature = 'gauss_newton' and yields a usable vcov", {
  res <- fit_fixture(keep = 1:8, curvature = "gauss_newton")
  fit <- res$fit
  expect_true(all(is.finite(coef(fit))))
  sm <- summary(fit)
  expect_true(all(is.finite(diag(sm$vcov))))
  expect_true(all(diag(sm$vcov) >= -1e-8))
})

test_that("curvature is validated against its allowed values", {
  cs <- curvature_surface()
  expect_error(
    terradish_algorithm(cs$f, leastsquares, cs$surface, cs$S, c(0, 0),
                        curvature = "nonsense"),
    "should be one of"
  )
})
