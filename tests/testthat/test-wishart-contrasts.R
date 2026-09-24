test_that("covariance likelihood and derivatives use site contrasts", {
  set.seed(261)
  n <- 6L
  Z <- matrix(rnorm(n * 12), n)
  E <- tcrossprod(Z) / 12
  S <- rWishart(1, 40, 0.8 * E + 0.2 * diag(n))[, , 1] / 40
  H <- diag(n) - matrix(1 / n, n, n)
  L <- qr.Q(qr(stats::contr.helmert(n)))
  phi <- c(tau = 0.8, sigma = log(0.2))
  fit <- wishart_covariance(E, S, phi, 40, validate = TRUE)
  centered <- wishart_covariance(E, H %*% S %*% H, phi, 40)
  shifted <- wishart_covariance(E, S + 3, phi, 40)
  distance <- outer(diag(S), diag(S), "+") - 2 * S
  generalized <- generalized_wishart(E, distance, phi, 40)
  expect_equal(fit$objective, centered$objective, tolerance = 1e-10)
  expect_equal(fit$objective, shifted$objective, tolerance = 1e-10)
  expect_equal(fit$objective, generalized$objective, tolerance = 1e-10)
  expect_equal(fit$gradient, generalized$gradient, tolerance = 1e-8)
  expect_equal(fit$hessian, generalized$hessian, tolerance = 1e-7)
  expect_equal(c(fit$gradient), fit$num_gradient, tolerance = 1e-6)
  expect_equal(unname(fit$hessian), fit$num_hessian, tolerance = 1e-5)
  expect_equal(fit$gradient_E, fit$num_gradient_E, tolerance = 1e-6)
  expect_equal(unname(fit$partial_E), fit$num_partial_E, tolerance = 1e-6)
  expect_equal(unname(fit$partial_S), fit$num_partial_S, tolerance = 1e-6)
  expect_equal(fit$jacobian_E(E), fit$num_jacobian_E(E), tolerance = 1e-6)
  expect_equal(fit$jacobian_S(E), fit$num_jacobian_S(E), tolerance = 1e-6)
  expect_equal(fit$fitted, 0.8 * E + 0.2 * diag(n))
  expect_equal(crossprod(L, fit$fitted %*% L),
               -0.5 * crossprod(L, generalized$fitted %*% L), tolerance = 1e-10)
})
