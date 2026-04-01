test_that("radish respects cores argument", {
  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)

  fit1 <- suppressWarnings(
    radish(melip.Fst ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
           cores = 1)
  )
  fit2 <- suppressWarnings(
    radish(melip.Fst ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE),
           cores = 2)
  )

  expect_equal(unname(coef(fit1)), unname(coef(fit2)), tolerance = 1e-6)
  expect_equal(unname(fit1$loglik), unname(fit2$loglik), tolerance = 1e-6)
})

test_that("radish_grid and radish_distance match across cores", {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("terradish"))
  {
    skip("parallel worker consistency is validated against the installed package")
  }

  dat <- melip_fixture(1:12)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta <- as.matrix(expand.grid(forestcover = c(-0.3, 0.3),
                                 altitude = c(-0.3, 0.3)))

  grid1 <- radish_grid(theta, melip.Fst ~ altitude + forestcover,
                       data = surface,
                       conductance_model = loglinear_conductance,
                       measurement_model = leastsquares,
                       cores = 1)
  grid2 <- radish_grid(theta, melip.Fst ~ altitude + forestcover,
                       data = surface,
                       conductance_model = loglinear_conductance,
                       measurement_model = leastsquares,
                       cores = 2)

  dist1 <- radish_distance(theta, ~ altitude + forestcover,
                           data = surface,
                           conductance_model = loglinear_conductance,
                           cores = 1)
  dist2 <- radish_distance(theta, ~ altitude + forestcover,
                           data = surface,
                           conductance_model = loglinear_conductance,
                           cores = 2)

  expect_equal(grid1$loglik, grid2$loglik, tolerance = 1e-8)
  expect_equal(dist1$distance, dist2$distance, tolerance = 1e-8)
})
