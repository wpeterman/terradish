# Two small contracts:
# 1. Measurement-model strings in terradish_cv_folds() resolve to the model
#    their name states, and "wishart" keeps its original meaning.
# 2. gaussian_scale_summary() warns about degrees once per session.

test_that("measurement-model strings resolve to the named model", {
  resolve <- terradish:::.terradish_measurement_model
  expect_identical(resolve("wishart_covariance"), wishart_covariance)
  expect_identical(resolve("generalized_wishart"), generalized_wishart)
  # Backward compatibility: the short string is unchanged.
  expect_identical(resolve("wishart"), generalized_wishart)
  expect_identical(resolve("mlpe"), mlpe)
  expect_identical(resolve("ls"), leastsquares)
  expect_error(resolve("covariance"), "Unknown measurement model")
  expect_error(resolve(1), "wishart_covariance")
})

test_that("the longitude/latitude sigma warning fires once per session", {
  skip_on_cran()
  dat <- melip_fixture(1:9)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates[["altitude"]], dat$coords,
                                 saveStack = TRUE)
  factory <- gaussian_smoothed_loglinear_conductance(surface)
  fit <- suppressWarnings(terradish(
    melip.Fst ~ altitude, surface, conductance_model = factory,
    measurement_model = mlpe, optimizer = "bfgs"))

  state <- terradish:::.terradish_notice_state
  old <- state$lonlat_sigma
  on.exit(state$lonlat_sigma <- old, add = TRUE)
  state$lonlat_sigma <- NULL

  expect_warning(first <- gaussian_scale_summary(fit), "once per session")
  expect_no_warning(second <- gaussian_scale_summary(fit))
  expect_identical(first, second)
  expect_identical(unique(first$native_unit), "degrees")

  # A supplied conversion never warns, whatever the session state.
  state$lonlat_sigma <- NULL
  expect_no_warning(gaussian_scale_summary(fit, distance_per_map_unit = 111,
                                           distance_unit = "km"))
})
