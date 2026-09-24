test_that("environmental kernels preserve pairwise differences and subset scales", {
  set.seed(911)
  x <- cbind(a = rnorm(9), b = rnorm(9))
  H <- diag(9) - matrix(1 / 9, 9, 9)
  for (transform in c("absdiff", "sqdiff", "euclidean", "manhattan")) {
    z <- pairwise_endpoint_covariates(x, transform = transform, scale = TRUE)
    g <- wishart_covariates(z)
    kernels <- attr(g, "kernel_covariates")
    for (k in seq_len(dim(kernels)[3])) {
      K <- kernels[, , k]
      expect_gte(min(eigen(K, symmetric = TRUE)$values), -1e-10)
      D <- outer(diag(K), diag(K), "+") - 2 * K
      expect_equal(c(D[lower.tri(D)]), c(z[, k]), tolerance = 1e-10)
    }
    for (index in list(c(2, 4, 6, 9), c(1, 2))) {
      subset <- attr(g, "subsetter")(index)
      direct <- wishart_covariates(attr(z, "site_covariates")[index, , drop = FALSE],
                                    transform = transform)
      expect_equal(attr(subset, "kernel_covariates"), attr(direct, "kernel_covariates"),
                    tolerance = 1e-10)
    }
    expect_s3_class(mlpe_covariates(z), "terradish_measurement_model")
  }
  squared <- attr(wishart_covariates(x, transform = "sqdiff"), "kernel_covariates")
  expect_equal(squared[, , 1], tcrossprod(c(H %*% x[, 1])), tolerance = 1e-10)
  bad <- structure(matrix(c(1, 1, 100), ncol = 1),
                    class = c("terradish_pairwise_covariates", "matrix", "array"))
  expect_error(wishart_covariates(bad), "not conditionally negative definite")
})

test_that("both environmental likelihood paths match dense site contrasts", {
  set.seed(912)
  n <- 8L
  x <- rnorm(n)
  E <- tcrossprod(matrix(rnorm(n * 12), n)) / 12
  gc <- wishart_covariates(x)
  gd <- wishart_covariates(x, model = "generalized_wishart")
  K <- attr(gc, "kernel_covariates")[, , 1]
  Sigma <- 0.8 * E + 0.6 * K + 0.2 * diag(n)
  S <- rWishart(1, 40, Sigma)[, , 1] / 40
  D <- outer(diag(S), diag(S), "+") - 2 * S
  phi <- c(0.8, 0.6, log(0.2))
  a <- gc(E, S, phi, 40)
  b <- gd(E, D, phi, 40)
  L <- qr.Q(qr(stats::contr.helmert(n)))
  V <- crossprod(L, Sigma %*% L)
  C <- crossprod(L, S %*% L)
  expected <- 20 * (as.numeric(determinant(V, logarithm = TRUE)$modulus) + sum(diag(solve(V, C))))
  expect_equal(a$objective, expected, tolerance = 1e-10)
  expect_equal(a$objective, b$objective, tolerance = 1e-10)
  expect_equal(a$gradient, b$gradient, tolerance = 1e-8)
  expect_equal(a$hessian, b$hessian, tolerance = 1e-8)
  expect_equal(c(a$gradient), numDeriv::grad(function(p) gc(E, S, p, 40)$objective, phi),
                tolerance = 1e-6)
  expect_equal(unname(a$hessian), numDeriv::hessian(function(p) gc(E, S, p, 40)$objective, phi),
                tolerance = 1e-5)
})
