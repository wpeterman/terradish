test_that("experimental Kron reduction matches the full reduced Laplacian block", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  conductance <- rep(1, nrow(surface$x))

  full_laplacian <- terradish:::.terradish_full_laplacian(surface, conductance)
  reduced_index <- terradish:::.graph_reduced_index(surface, length(conductance))
  reduced_laplacian <- terradish:::.graph_reduced_laplacian(surface, conductance)

  laplacian_diff <- full_laplacian[reduced_index, reduced_index] - reduced_laplacian
  max_abs_diff <- if (length(laplacian_diff@x)) max(abs(laplacian_diff@x)) else 0
  expect_equal(max_abs_diff, 0, tolerance = 1e-10)
})

test_that("experimental Kron reduction preserves focal effective distances", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  theta <- c(-0.3, 0.3)
  conductance <- model(theta)$conductance

  full_fit <- terradish_algorithm(
    model,
    leastsquares,
    surface,
    diag(length(surface$demes)),
    nu = 1000,
    theta = theta,
    objective = FALSE,
    gradient = FALSE,
    hessian = FALSE,
    partial = FALSE,
    solver = "direct"
  )

  kron <- terradish_kron_reduce(surface, conductance, covariance = TRUE)

  expect_s3_class(kron, "terradish_kron_reduction")
  expect_equal(kron$n_boundary, length(surface$demes))
  expect_equal(
    dist_from_cov(kron$covariance),
    dist_from_cov(as.matrix(full_fit$covariance)),
    tolerance = 1e-6
  )
})
