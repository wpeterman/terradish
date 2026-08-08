test_that("check_distance_response identifies squared-Euclidean responses", {
  xy <- matrix(c(0, 0,
                 1, 0,
                 0, 1,
                 1, 1), ncol = 2, byrow = TRUE)
  D <- as.matrix(dist(xy))^2

  checked <- check_distance_response(D)

  expect_named(
    checked,
    c("admissible", "min_eigenvalue", "relative_min_eigenvalue",
      "n_negative", "tolerance", "eigenvalue_threshold",
      "max_eigenvalue", "n_sites")
  )
  expect_true(checked$admissible)
  expect_equal(checked$n_negative, 0L)
  expect_equal(checked$n_sites, 4L)
  expect_equal(checked$tolerance, 1e-7)
})

test_that("check_distance_response tolerates numerical-scale violations", {
  xy <- matrix(c(0, 0,
                 1, 0,
                 0, 1,
                 1, 1), ncol = 2, byrow = TRUE)
  D <- as.matrix(dist(xy))^2
  D[1, 2] <- D[2, 1] <- D[1, 2] - 1e-10

  checked <- check_distance_response(D, tol = 1e-7)

  expect_lt(checked$min_eigenvalue, 0)
  expect_true(checked$admissible)
  expect_equal(checked$n_negative, 0L)
})

test_that("check_distance_response reports substantive non-Euclidean geometry", {
  D <- matrix(c(0, 1, 1,
                1, 0, 9,
                1, 9, 0), nrow = 3, byrow = TRUE)

  checked <- check_distance_response(D)

  expect_false(checked$admissible)
  expect_lt(checked$min_eigenvalue, 0)
  expect_lt(checked$relative_min_eigenvalue, -1e-7)
  expect_gt(checked$n_negative, 0L)
})

test_that("check_distance_response rejects malformed responses", {
  expect_error(check_distance_response(1:4), "numeric square matrix")
  expect_error(check_distance_response(matrix(letters[1:4], 2)), "numeric")
  expect_error(check_distance_response(matrix(1:6, 2, 3)), "square")
  expect_error(check_distance_response(matrix(0, 1, 1)), "at least two")

  nonfinite <- matrix(c(0, Inf, Inf, 0), 2)
  expect_error(check_distance_response(nonfinite), "finite")

  asymmetric <- matrix(c(0, 1, 2, 0), 2)
  expect_error(check_distance_response(asymmetric), "symmetric")

  nonzero_diagonal <- matrix(c(0.1, 1, 1, 0), 2)
  expect_error(check_distance_response(nonzero_diagonal), "zero diagonal")

  negative <- matrix(c(0, -1, -1, 0), 2)
  expect_error(check_distance_response(negative), "nonnegative")

  expect_error(check_distance_response(diag(2), tol = 0), "positive")
})

test_that("generalized-Wishart fitting rejects an inadmissible response", {
  E <- diag(3)
  D <- matrix(c(0, 1, 1,
                1, 0, 9,
                1, 9, 0), nrow = 3, byrow = TRUE)

  expect_error(
    generalized_wishart(E, D, phi = c(tau = 1, sigma = 0), nu = 20),
    "not squared Euclidean"
  )

  dat <- melip_fixture(keep = 1:4)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  D4 <- as.matrix(dist(matrix(c(0, 0,
                                1, 0,
                                0, 1,
                                1, 1), ncol = 2, byrow = TRUE)))^2
  D4[1, 2] <- D4[2, 1] <- 100

  expect_error(
    terradish(D4 ~ altitude + forestcover,
              data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = generalized_wishart,
              nu = 20),
    "not squared Euclidean"
  )
})
