# Does curvature = "gauss_newton" give the same standard errors as the default
# exact curvature for an MLPE fit?
#
# Short answer recorded by these tests: the coefficients agree, the standard
# errors do not in general. Gauss-Newton drops the residual-weighted second
# derivative terms, and an MLPE fit to real FST data has residuals that are
# not negligible. On the melip data the two sets of standard errors differ by
# a few percent on a small subset and by roughly a factor of two on the full
# set of sites. The numbers in the comments are the values observed when the
# tests were written (R 4.6.2, Windows); the tolerances were chosen from
# them and were not adjusted to make a comparison pass.

mlpe_curvature_fixture <- function(keep, aggregate_factor) {
  dat <- melip_fixture(keep, aggregate_factor = aggregate_factor)
  # Keep one site per raster cell so the graph has no co-located sites.
  one_per_cell <- !duplicated(terra::cellFromXY(dat$covariates,
                                                terra::crds(dat$coords)))
  list(Fst = dat$melip.Fst[one_per_cell, one_per_cell],
       surface = conductance_surface(dat$covariates,
                                     dat$coords[one_per_cell]))
}

standard_errors <- function(hessian) sqrt(diag(solve(hessian)))

test_that("gauss_newton is supported for mlpe and reproduces the exact coefficients", {
  fx <- mlpe_curvature_fixture(keep = 1:12, aggregate_factor = 3)
  Fst <- fx$Fst

  exact <- terradish(Fst ~ altitude + forestcover, fx$surface,
                     measurement_model = mlpe, curvature = "exact")
  gauss_newton <- terradish(Fst ~ altitude + forestcover, fx$surface,
                            measurement_model = mlpe,
                            curvature = "gauss_newton")
  expect_equal(exact$convergence$code, 0)
  expect_equal(gauss_newton$convergence$code, 0)

  # Same optimum. The optimizer takes different steps under the two
  # curvatures, so agreement is to optimizer tolerance, not bit for bit.
  # Observed largest absolute coefficient difference: 6.8e-07.
  # Observed log-likelihood difference: 5.4e-11.
  expect_equal(coef(gauss_newton), coef(exact), tolerance = 1e-4)
  expect_equal(gauss_newton$loglik, exact$loglik, tolerance = 1e-8)

  # Standard errors: observed exact 0.0894 (altitude) and 0.0989
  # (forestcover); Gauss-Newton 0.0932 and 0.1060. Relative differences
  # +4.3% and +7.2%. The bound is about twice the larger observed value.
  se_exact <- sqrt(diag(vcov(exact)))
  se_gauss_newton <- sqrt(diag(vcov(gauss_newton)))
  expect_true(all(is.finite(se_gauss_newton)) && all(se_gauss_newton > 0))
  relative_difference <- se_gauss_newton / se_exact - 1
  expect_lt(max(abs(relative_difference)), 0.15)
  # They are close here but they are not the same quantity.
  expect_gt(max(abs(relative_difference)), 0.01)
})

test_that("at one optimum the two curvatures share the likelihood but not the standard errors", {
  skip_on_cran()
  # All melip sites that occupy distinct cells of the coarsened raster (31).
  fx <- mlpe_curvature_fixture(keep = NULL, aggregate_factor = 3)
  Fst <- fx$Fst

  exact <- terradish(Fst ~ altitude + forestcover, fx$surface,
                     measurement_model = mlpe, curvature = "exact")
  expect_equal(exact$convergence$code, 0)

  # Evaluate both curvatures at the same parameter values, so the comparison
  # does not depend on the path an optimizer took.
  model <- loglinear_conductance(~ altitude + forestcover, fx$surface$x)
  evaluate <- function(curvature)
    terradish_algorithm(model, mlpe, fx$surface, Fst, theta = coef(exact),
                        partial = FALSE, curvature = curvature)
  at_exact <- evaluate("exact")
  at_gauss_newton <- evaluate("gauss_newton")

  # The curvature option changes only the Hessian.
  expect_equal(at_gauss_newton$objective, at_exact$objective)
  expect_equal(at_gauss_newton$gradient, at_exact$gradient)

  # Both are positive definite, so both give usable standard errors.
  expect_true(all(eigen(at_exact$hessian, symmetric = TRUE,
                        only.values = TRUE)$values > 0))
  expect_true(all(eigen(at_gauss_newton$hessian, symmetric = TRUE,
                        only.values = TRUE)$values > 0))

  # The exact Hessian here is the one behind vcov() of the fit.
  se_exact <- standard_errors(at_exact$hessian)
  expect_equal(unname(se_exact), unname(sqrt(diag(vcov(exact)))),
               tolerance = 1e-6)

  # Observed: exact 0.150 (altitude) and 0.193 (forestcover); Gauss-Newton
  # 0.303 and 0.341. Relative differences +102% and +77%. This is a material
  # disagreement, so the test records it instead of asserting agreement:
  # Gauss-Newton standard errors must not be reported as interchangeable
  # with the exact ones for MLPE fits. If a later change brings the two
  # within 25% on this example, revisit the documentation of `curvature`.
  se_gauss_newton <- standard_errors(at_gauss_newton$hessian)
  relative_difference <- se_gauss_newton / se_exact - 1
  expect_gt(max(abs(relative_difference)), 0.25)
})
