test_that("landmark approximation subsets focal populations and rhs columns", {
  dat <- melip_fixture(1:10)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  subset <- terradish:::.terradish_landmark_subset(
    surface,
    dat$melip.Fst,
    approximation_control = list(
      n_landmarks = 4L,
      method = "spacefill",
      seed = 1
    )
  )

  expect_true(subset$used)
  expect_length(subset$index, 4L)
  expect_equal(nrow(subset$S), 4L)
  expect_equal(ncol(subset$S), 4L)
  expect_equal(length(subset$data$demes), 4L)
  expect_equal(ncol(subset$data$rhs), 4L)
})

test_that("terradish landmark approximation refines back to the full-data fit for mlpe", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  control <- NewtonRaphsonControl(maxit = 4, verbose = FALSE)

  fit_exact <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct"
    )
  )

  fit_landmark <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = mlpe,
      control = control,
      solver = "direct",
      approximation = "landmark",
      approximation_control = list(
        n_landmarks = 4L,
        method = "spacefill",
        seed = 1,
        exact_refine = TRUE
      )
    )
  )

  expect_true(isTRUE(fit_landmark$approximation$used))
  expect_equal(fit_landmark$approximation$type, "landmark")
  expect_equal(fit_landmark$approximation$n_landmarks, 4L)
  expect_equal(dim(fit_landmark$fit$covariance), dim(fit_exact$fit$covariance))
  expect_equal(fit_landmark$loglik, fit_exact$loglik, tolerance = 1e-4)
  expect_lt(max(abs(fit_landmark$mle$theta - fit_exact$mle$theta)), 2e-2)
})

test_that("landmark exact refinement is guarded for leastsquares", {
  dat <- melip_fixture(1:8)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  control <- NewtonRaphsonControl(maxit = 4, verbose = FALSE)

  fit_landmark <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = leastsquares,
      control = control,
      solver = "direct",
      approximation = "landmark",
      approximation_control = list(
        n_landmarks = 4L,
        method = "spacefill",
        seed = 1,
        exact_refine = TRUE
      )
    )
  )

  expect_equal(fit_landmark$approximation$type, "landmark")
  expect_false(isTRUE(fit_landmark$approximation$used))
  expect_equal(fit_landmark$approximation$refine_guard, "disabled_for_leastsquares")
})
