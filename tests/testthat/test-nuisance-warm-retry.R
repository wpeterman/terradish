test_that("a finite failed nuisance warm start retries the default start", {
  f <- readRDS(test_path("fixtures", "wishart-warm-nonconvergence.rds"))
  g <- wishart_covariates(f$pairs, model = "wishart_covariance")
  expect_true(f$original_code != 0L)
  recovered <- expect_fixture_warnings(
    radish_subproblem(g, f$E, f$S, f$nu, phi = f$phi, control = f$control),
    "`maxit` reached for Newton steps")
  cold <- radish_subproblem(g, f$E, f$S, f$nu, control = f$control)
  expect_equal(recovered$convergence, 0L)
  expect_equal(cold$convergence, 0L)
  expect_equal(recovered$loglikelihood, cold$loglikelihood, tolerance = 1e-8)
  expect_equal(recovered$phi, cold$phi, tolerance = 1e-7)
  expect_true(all(c(recovered$phi) >= f$phi$lower))
  expect_true(all(c(recovered$phi) <= f$phi$upper))
})

test_that("an unsuccessful default retry does not clear nuisance failure", {
  f <- readRDS(test_path("fixtures", "wishart-warm-nonconvergence.rds"))
  g <- wishart_covariates(f$pairs, model = "wishart_covariance")
  f$control$maxit <- 1L
  fit <- expect_fixture_warnings(
    radish_subproblem(g, f$E, f$S, f$nu, phi = f$phi, control = f$control),
    "`maxit` reached for Newton steps")
  expect_true(fit$convergence != 0L)
})
