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

test_that("tiled Kron reduction equals the single-shot reduction exactly", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance

  single <- terradish_kron_reduce(surface, conductance)
  tiled  <- terradish_kron_reduce_tiled(surface, conductance, n_tiles = 16L)

  expect_s3_class(tiled, "terradish_kron_reduction")
  diff <- as.matrix(single$laplacian) - as.matrix(tiled$laplacian[
    match(single$boundary, tiled$boundary),
    match(single$boundary, tiled$boundary)])
  expect_lt(max(abs(diff)) / max(abs(as.matrix(single$laplacian))), 1e-8)
})

test_that("tiled Kron reduction is exact regardless of tile count", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.2, 0.4))$conductance

  a <- terradish_kron_reduce_tiled(surface, conductance, n_tiles = 4L)
  b <- terradish_kron_reduce_tiled(surface, conductance, n_tiles = 25L)
  expect_equal(as.matrix(a$laplacian), as.matrix(b$laplacian), tolerance = 1e-8)
  # more tiles means a smaller largest interior solve (lower peak memory)
  expect_lte(b$peak$interior, a$peak$interior)
})

test_that("tiled reduction preserves focal effective distances and bounds memory", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  theta <- c(-0.3, 0.3)
  conductance <- model(theta)$conductance

  full_fit <- terradish_algorithm(
    model, leastsquares, surface, diag(length(surface$demes)),
    nu = 1000, theta = theta, objective = FALSE, gradient = FALSE,
    hessian = FALSE, partial = FALSE, solver = "direct")

  tiled <- terradish_kron_reduce_tiled(surface, conductance, covariance = TRUE)
  expect_equal(
    dist_from_cov(tiled$covariance),
    dist_from_cov(as.matrix(full_fit$covariance)),
    tolerance = 1e-6)
  # the whole point: the largest interior solve is far smaller than eliminating
  # the entire interior at once
  expect_lt(tiled$peak$interior, tiled$n_interior)
})

test_that("tiled reduction accepts a user-supplied vertex partition", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(0, 0))$conductance

  n <- nrow(surface$x)
  labels <- ((seq_len(n) - 1L) %% 6L) + 1L          # 6 interleaved tiles
  user   <- terradish_kron_reduce_tiled(surface, conductance, tiles = labels)
  single <- terradish_kron_reduce(surface, conductance)
  expect_equal(
    as.matrix(user$laplacian[match(single$boundary, user$boundary),
                             match(single$boundary, user$boundary)]),
    as.matrix(single$laplacian),
    tolerance = 1e-8)

  expect_error(
    terradish_kron_reduce_tiled(surface, conductance, tiles = labels[-1]),
    "one entry per graph vertex")
})

test_that("parallel tiled reduction matches the sequential and single-shot results", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance

  seq1 <- terradish_kron_reduce_tiled(surface, conductance, n_tiles = 16L, cores = 1L)
  par2 <- terradish_kron_reduce_tiled(surface, conductance, n_tiles = 16L, cores = 2L)
  single <- terradish_kron_reduce(surface, conductance)

  # parallel == sequential, exactly
  expect_equal(as.matrix(par2$laplacian), as.matrix(seq1$laplacian), tolerance = 1e-10)
  expect_equal(par2$cores, 2)
  # and both match the single-shot reduction
  P <- match(single$boundary, seq1$boundary)
  expect_lt(
    max(abs(as.matrix(single$laplacian) - as.matrix(seq1$laplacian[P, P]))) /
      max(abs(as.matrix(single$laplacian))),
    1e-8)
})
