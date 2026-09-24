# Establish F2's parameterization before changing the covariance likelihood.
pkgload::load_all()
set.seed(20260924)
for (n in c(4L, 9L)) {
  X <- matrix(rnorm(n * n), n)
  E <- tcrossprod(X) / n + diag(n)
  Sigma <- 0.7 * E + 0.2 * diag(n)
  C <- rWishart(1, 30, Sigma)[, , 1] / 30
  L <- qr.Q(qr(stats::contr.helmert(n)))
  contrast_sigma <- crossprod(L, Sigma %*% L)
  contrast_response <- crossprod(L, C %*% L)
  objective <- 30 / 2 * (as.numeric(determinant(contrast_sigma)$modulus) +
    sum(diag(solve(contrast_sigma, contrast_response))))
  phi <- c(tau = 0.7, sigma = log(0.2))
  gw <- generalized_wishart(E, dist_from_cov(C), phi = phi, nu = 30)
  expected_distance <- 0.7 * dist_from_cov(E) +
    0.4 * (matrix(1, n, n) - diag(n))
  cat("n =", n, "objective difference =", objective - gw$objective,
      "distance difference =", max(abs(expected_distance - gw$fitted)), "\n")
  stopifnot(abs(objective - gw$objective) < 1e-10,
            max(abs(expected_distance - gw$fitted)) < 1e-10)
}
cat("Same tau and log(sigma) in both parameterizations; off-diagonal nugget is 2*sigma.\n")
