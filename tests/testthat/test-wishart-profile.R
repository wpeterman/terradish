test_that("a vanishing warm nugget cannot retain an inferior Wishart mode", {
  b <- readRDS(test_path("fixtures", "wishart-nuisance-mode.rds"))
  profile <- function(phi = NULL, S = b$S, g = wishart_covariance)
    radish_subproblem(g, b$E, S, nu = b$nu, phi = phi,
      control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10))
  cold <- profile()
  warm <- profile(b$retained_phi)
  expect_equal(cold$convergence, 0)
  expect_equal(warm$convergence, 0)
  expect_equal(warm$loglikelihood, cold$loglikelihood, tolerance = 1e-12)
  expect_equal(warm$phi, cold$phi, tolerance = 1e-12)
  expect_equal(-cold$loglikelihood, -151494.469917915, tolerance = 1e-10)
  expect_gt(-cold$loglikelihood - (-151595.12989754), 100)
  expect_gt(exp(warm$phi[2]), .1)
  expect_lt(max(abs(warm$fit$gradient)), 1e-6)
  # Evaluation history cannot change a repeated objective or its derivatives.
  other <- profile(c(1, -80))
  expect_equal(other$gradient, cold$gradient, tolerance = 1e-12)
  expect_equal(other$fit$hessian, cold$fit$hessian, tolerance = 1e-12)
  distance <- profile(S = dist_from_cov(b$S), g = generalized_wishart)
  expect_equal(unname(distance$phi), unname(cold$phi), tolerance = 1e-6)
  expect_equal(distance$loglikelihood, cold$loglikelihood, tolerance = 1e-10)
})

test_that("profile derivatives agree with perturbations of the observed case", {
  b <- readRDS(test_path("fixtures", "wishart-nuisance-mode.rds"))
  evaluate <- function(E) radish_subproblem(wishart_covariance, E, b$S,
    nu = b$nu, phi = b$retained_phi,
    control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10))
  D <- diag(seq(.001, .002, length.out = nrow(b$E)))
  fit <- evaluate(b$E)
  step <- 1e-3
  plus <- evaluate(b$E + step * D)
  minus <- evaluate(b$E - step * D)
  numerical_gradient <- (plus$loglikelihood - minus$loglikelihood) / (2 * step)
  expect_equal(sum(fit$gradient * D), numerical_gradient, tolerance = 1e-5)
  numerical_hessian <- (plus$gradient - minus$gradient) / (2 * step)
  expect_equal(fit$jacobian_E(D), numerical_hessian, tolerance = 1e-5)
})

test_that("Wishart profiling preserves boundary and fixed-parameter contracts", {
  E <- diag(c(.1, .4, 1, 3))
  pure_nugget <- radish_subproblem(wishart_covariance, E, diag(4), nu = 30)
  expect_equal(pure_nugget$phi[1, 1], 0, ignore_attr = TRUE)
  expect_true(pure_nugget$boundary)
  pure_kernel <- radish_subproblem(wishart_covariance, E, .8 * E, nu = 30)
  expect_equal(pure_kernel$phi[1, 1], .8, tolerance = 1e-7, ignore_attr = TRUE)
  expect_lt(exp(pure_kernel$phi[2, 1]), 1e-12)
  fixed <- c(tau = .3, sigma = log(.4))
  fit <- radish_subproblem(wishart_covariance, E, diag(4), nu = 30,
    phi = list(phi = fixed, lower = fixed, upper = fixed))
  expect_equal(c(fit$phi), unname(fixed))
  expect_equal(fit$loglikelihood,
    wishart_covariance(E, diag(4), fixed, nu = 30)$objective)
  bounded <- radish_subproblem(wishart_covariance, E, diag(4), nu = 30,
    phi = list(phi = c(.4, 0), lower = c(.4, -Inf), upper = c(1, Inf)))
  expect_gte(bounded$phi[1, 1], .4)
})
