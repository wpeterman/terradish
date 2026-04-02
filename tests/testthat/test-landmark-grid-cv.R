test_that("terradish_grid records landmark approximation metadata", {
  fx <- fit_fixture(keep = 1:8,
                    control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst
  theta0 <- coef(fx$fit)
  theta <- rbind(theta0, theta0 + rep(c(0.05, -0.05), length.out = length(theta0)))
  colnames(theta) <- names(theta0)

  grid <- suppressWarnings(
    terradish_grid(
      theta = theta,
      formula = melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      cores = 1,
      approximation = "landmark",
      approximation_control = list(
        n_landmarks = 4L,
        method = "spacefill",
        seed = 1
      )
    )
  )

  expect_s3_class(grid, "terradish_grid")
  expect_s3_class(grid, "radish_grid")
  expect_equal(grid$approximation$type, "landmark")
  expect_true(isTRUE(grid$approximation$used))
  expect_equal(grid$approximation$n_landmarks, 4L)
  expect_length(grid$loglik, nrow(theta))
})

test_that("terradish_grid coarse-raster screening keeps the full focal set", {
  fx <- fit_fixture(keep = 1:12,
                    control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst
  theta0 <- coef(fx$fit)
  theta <- rbind(theta0, theta0 + rep(c(0.05, -0.05), length.out = length(theta0)))
  colnames(theta) <- names(theta0)

  grid <- suppressWarnings(
    terradish_grid(
      theta = theta,
      formula = melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      cores = 1,
      approximation = "coarse_raster",
      approximation_control = list(factor = 2L)
    )
  )

  expect_s3_class(grid, "terradish_grid")
  expect_s3_class(grid, "radish_grid")
  expect_equal(grid$approximation$type, "coarse_raster")
  expect_true(isTRUE(grid$approximation$used))
  expect_equal(grid$approximation$full_focal, nrow(melip.Fst))
  expect_lt(grid$approximation$coarse_vertices, grid$approximation$full_vertices)
  expect_length(grid$loglik, nrow(theta))
  expect_true(all(is.finite(grid$loglik)))
})

test_that("terradish_cv forwards landmark approximation to held-out grid evaluation", {
  dat <- melip_fixture(1:16)
  melip.Fst <- dat$melip.Fst
  approximation_control <- list(
    n_landmarks = 4L,
    method = "spacefill",
    seed = 11
  )

  out <- suppressWarnings(
    terradish_cv(
      dat$coords,
      dat$covariates,
      melip.Fst ~ altitude + forestcover,
      model = "mlpe",
      prop_train = 0.70,
      seed = 3,
      fit_full = FALSE,
      control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
      approximation = "landmark",
      approximation_control = approximation_control
    )
  )

  expect_true(isTRUE(out$train_mod$approximation$used))
  expect_equal(out$train_mod$approximation$type, "landmark")

  melip.Fst_test <- dat$melip.Fst[out$test_index, out$test_index, drop = FALSE]
  test_surface <- conductance_surface(
    dat$covariates,
    dat$coords[out$test_index],
    directions = 8
  )

  manual <- suppressWarnings(
    terradish_grid(
      theta = matrix(coef(out$train_mod), nrow = 1),
      formula = melip.Fst_test ~ altitude + forestcover,
      data = test_surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      cores = 1,
      approximation = "landmark",
      approximation_control = approximation_control
    )
  )

  expect_true(isTRUE(manual$approximation$used))
  expect_equal(out$cv_loglik, manual$loglik)
})

test_that("terradish_cv can use coarse-raster screening on the held-out grid only", {
  dat <- melip_fixture(1:16)
  melip.Fst <- dat$melip.Fst

  out <- suppressWarnings(
    terradish_cv(
      dat$coords,
      dat$covariates,
      melip.Fst ~ altitude + forestcover,
      model = "ls",
      prop_train = 0.70,
      seed = 3,
      fit_full = FALSE,
      control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
      approximation = "coarse_raster",
      approximation_control = list(factor = 2L)
    )
  )

  expect_true(is.list(out))
  expect_length(out$cv_loglik, 1L)
  expect_equal(out$train_mod$approximation$type, "none")
  expect_false(isTRUE(out$train_mod$approximation$used))
})
