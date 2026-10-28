test_that("grid and distance evaluation preserve the caller's session settings", {
  dat <- melip_fixture(1:6, aggregate_factor = 5L)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  theta <- matrix(c(-0.2, 0.2), nrow = 1,
                  dimnames = list(NULL, c("altitude", "forestcover")))
  response <- dat$melip.Fst
  wd <- getwd()
  opts <- options()

  grid <- terradish_grid(theta, response ~ altitude + forestcover, surface,
                         loglinear_conductance, leastsquares, cores = 1)
  distance <- terradish_distance(theta, ~ altitude + forestcover, surface,
                                 loglinear_conductance, cores = 1)

  expect_identical(getwd(), wd)
  expect_identical(options(), opts)
  expect_true(all(is.finite(grid$loglik)))
  expect_true(all(is.finite(distance$distance)))

  expect_error(terradish_distance(theta[, 1, drop = FALSE],
                                  ~ altitude + forestcover, surface,
                                  loglinear_conductance, cores = 1))
  expect_identical(getwd(), wd)
  expect_identical(options(), opts)
})
