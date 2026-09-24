test_that("MLPE correlation algebra agrees with dense matrices near both bounds", {
  set.seed(410)
  pairs <- which(lower.tri(diag(6)), arr.ind = TRUE)
  U <- Matrix::sparseMatrix(i = rep(seq_len(nrow(pairs)), 2), j = c(pairs), x = 1)
  x <- matrix(rnorm(nrow(U) * 2), ncol = 2)
  for (logit in c(-20, -2, 2, 10)) {
    rho <- .5 * plogis(logit)
    Sigma <- (1 - 2 * rho) * diag(nrow(U)) + rho * as.matrix(tcrossprod(U))
    operator <- .mlpe_correlation_operator(U, logit)
    expect_equal(operator$inverse(x), solve(Sigma, x), tolerance = 1e-7)
    expect_equal(operator$quadratic(x), sum(x * solve(Sigma, x)), tolerance = 1e-7)
    expect_equal(operator$logdet, as.numeric(determinant(Sigma, logarithm = TRUE)$modulus), tolerance = 1e-7)
  }
  operator <- .mlpe_correlation_operator(U, 60)
  expect_gt(operator$quadratic(x), 0)
  expect_true(is.finite(operator$quadratic(x)))
})

test_that("an extreme upper-correlation trial cannot gain spurious likelihood", {
  set.seed(411)
  E <- crossprod(matrix(rnorm(36), 6))
  S <- as.matrix(dist(matrix(rnorm(12), 6)))
  base <- mlpe(E, S, phi = c(0, 1, 0, 60))
  extended <- mlpe_covariates(seq_len(6))(E, S, phi = c(0, 1, 0, 0, 60))
  expect_gt(base$objective, 1e12)
  expect_gt(extended$objective, 1e12)
})
