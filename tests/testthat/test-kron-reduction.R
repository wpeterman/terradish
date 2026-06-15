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

test_that("auto default uses the single-shot reduction and matches it", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance

  single <- terradish_kron_reduce(surface, conductance)
  tiled  <- terradish_kron_reduce_tiled(surface, conductance)   # auto -> direct (small graph)

  expect_s3_class(tiled, "terradish_kron_reduction")
  expect_equal(tiled$method, "direct")
  diff <- as.matrix(single$laplacian) - as.matrix(tiled$laplacian[
    match(single$boundary, tiled$boundary),
    match(single$boundary, tiled$boundary)])
  expect_lt(max(abs(diff)) / max(abs(as.matrix(single$laplacian))), 1e-8)
})

test_that("auto switches to nested dissection under a tight memory budget", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance
  single <- terradish_kron_reduce(surface, conductance)

  tiled <- terradish_kron_reduce_tiled(surface, conductance, mem_budget = 1)
  expect_equal(tiled$method, "nested")
  P <- match(single$boundary, tiled$boundary)
  expect_lt(
    max(abs(as.matrix(single$laplacian) - as.matrix(tiled$laplacian[P, P]))) /
      max(abs(as.matrix(single$laplacian))), 1e-8)
})

test_that("nested dissection is exact for any leaf size and bounds the factorization", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.2, 0.4))$conductance

  small <- terradish_kron_reduce_tiled(surface, conductance, method = "nested", n_tiles = 1000L)
  large <- terradish_kron_reduce_tiled(surface, conductance, method = "nested", n_tiles = 4000L)
  expect_equal(as.matrix(small$laplacian), as.matrix(large$laplacian), tolerance = 1e-8)
  # the largest single factorization is bounded well below the whole interior,
  # and a smaller leaf does not enlarge it
  expect_lt(small$peak$interior, small$n_interior)
  expect_lte(small$peak$interior, large$peak$interior)
})

test_that("nested dissection preserves focal effective distances", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  theta <- c(-0.3, 0.3)
  conductance <- model(theta)$conductance

  full_fit <- terradish_algorithm(
    model, leastsquares, surface, diag(length(surface$demes)),
    nu = 1000, theta = theta, objective = FALSE, gradient = FALSE,
    hessian = FALSE, partial = FALSE, solver = "direct")

  tiled <- terradish_kron_reduce_tiled(surface, conductance, method = "nested", covariance = TRUE)
  expect_equal(
    dist_from_cov(tiled$covariance),
    dist_from_cov(as.matrix(full_fit$covariance)),
    tolerance = 1e-6)
  expect_lt(tiled$peak$interior, tiled$n_interior)
})

test_that("explicit-partition path matches single-shot, parallel equals sequential", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance
  single <- terradish_kron_reduce(surface, conductance)

  n <- nrow(surface$x)
  labels <- ((seq_len(n) - 1L) %% 6L) + 1L          # 6 interleaved tiles
  seq1 <- terradish_kron_reduce_tiled(surface, conductance, tiles = labels, cores = 1L)
  par2 <- terradish_kron_reduce_tiled(surface, conductance, tiles = labels, cores = 2L)

  P <- match(single$boundary, seq1$boundary)
  expect_lt(
    max(abs(as.matrix(single$laplacian) - as.matrix(seq1$laplacian[P, P]))) /
      max(abs(as.matrix(single$laplacian))), 1e-8)
  expect_equal(as.matrix(par2$laplacian), as.matrix(seq1$laplacian), tolerance = 1e-10)
  expect_equal(par2$cores, 2)

  expect_error(
    terradish_kron_reduce_tiled(surface, conductance, tiles = labels[-1]),
    "one entry per graph vertex")
})

test_that("nested-dissection result is identical with cores > 1", {
  dat <- melip_fixture(1:8)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  model <- loglinear_conductance(~ altitude + forestcover, surface$x)
  conductance <- model(c(-0.3, 0.3))$conductance

  one <- terradish_kron_reduce_tiled(surface, conductance, method = "nested", cores = 1L)
  par <- terradish_kron_reduce_tiled(surface, conductance, method = "nested", cores = 2L)
  expect_equal(as.matrix(par$laplacian), as.matrix(one$laplacian), tolerance = 1e-10)
  expect_equal(par$cores, 2)
})
