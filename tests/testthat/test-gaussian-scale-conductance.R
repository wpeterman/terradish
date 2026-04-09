test_that("gaussian scale-aware conductance matches numerical derivatives for interactions and polynomials", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)

  factory <- gaussian_smoothed_loglinear_conductance(
    surface,
    sigma_lower = 0.2,
    sigma_upper = 5
  )
  conductance_model <- factory(~ altitude * forestcover + I(altitude^2), surface$x)
  theta <- attr(conductance_model, "default")
  theta[] <- c(
    altitude = 0.15,
    forestcover = -0.10,
    `altitude:forestcover` = 0.04,
    `I(altitude^2)` = 0.03,
    `sigma.altitude` = 1.2,
    `sigma.forestcover` = 0.8
  )[names(theta)]

  fit <- terradish_algorithm(
    conductance_model,
    leastsquares,
    surface,
    melip.Fst,
    theta = theta,
    gradient = TRUE,
    hessian = TRUE,
    partial = FALSE
  )

  num_gradient <- terradish:::.numderiv_grad(
    function(par) terradish_algorithm(
      conductance_model,
      leastsquares,
      surface,
      melip.Fst,
      theta = par,
      gradient = FALSE,
      hessian = FALSE,
      partial = FALSE
    )$objective,
    theta
  )
  num_hessian <- terradish:::.numderiv_hessian(
    function(par) terradish_algorithm(
      conductance_model,
      leastsquares,
      surface,
      melip.Fst,
      theta = par,
      gradient = FALSE,
      hessian = FALSE,
      partial = FALSE
    )$objective,
    theta
  )

  expect_equal(c(fit$gradient), c(num_gradient), tolerance = 1e-4)
  expect_equal(unname(fit$hessian), unname(num_hessian), tolerance = 2e-3)
})

test_that("gaussian scale-aware conductance fits jointly with bounded sigma", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["altitude"]]
  names(covariates) <- "altitude"
  surface <- conductance_surface(covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)

  factory <- gaussian_smoothed_loglinear_conductance(
    surface,
    sigma_lower = 0.5,
    sigma_upper = 2
  )
  conductance_model <- factory(~ altitude + I(altitude^2), surface$x)
  expect_equal(attr(conductance_model, "lower")[["sigma.altitude"]], 0.5)
  expect_equal(attr(conductance_model, "upper")[["sigma.altitude"]], 2)

  warnings <- character()
  fit <- withCallingHandlers(
    terradish(
      melip.Fst ~ altitude + I(altitude^2),
      data = surface,
      conductance_model = factory,
      measurement_model = leastsquares,
      optimizer = "auto",
      leverage = TRUE,
      control = NewtonRaphsonControl(maxit = 8, verbose = FALSE)
    ),
    warning = function(w)
    {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("Leverage diagnostics are not available", warnings)))

  expect_s3_class(fit, "terradish")
  expect_true(all(c("altitude", "I(altitude^2)", "sigma.altitude") %in%
                    names(coef(fit))))
  expect_gte(unname(coef(fit)[["sigma.altitude"]]), 0.5)
  expect_lte(unname(coef(fit)[["sigma.altitude"]]), 2)
  expect_null(fit$leverage$S)

  cond <- conductance(surface, fit)
  expect_s4_class(cond, "SpatRaster")

  smry <- suppressWarnings(summary(fit))
  expect_true(all(c("altitude", "I(altitude^2)", "sigma.altitude") %in%
                    rownames(smry$ztable)))

  sigma_table <- gaussian_scale_summary(
    fit,
    probabilities = c(0.5, 0.95),
    distance_per_map_unit = 0.001,
    distance_unit = "km"
  )
  expect_true(all(c("covariate", "sigma", "sigma_cells_x", "radial_95",
                    "sigma_km", "radial_95_km") %in% names(sigma_table)))
  expect_equal(sigma_table$covariate, "altitude")
})

test_that("gaussian scale-aware conductance validates sigma starts and keeps coarse raster disabled", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["altitude"]]
  names(covariates) <- "altitude"
  surface <- conductance_surface(covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)

  factory <- gaussian_smoothed_loglinear_conductance(
    surface,
    sigma_lower = 0.5,
    sigma_upper = 2
  )

  expect_error(
    terradish(
      melip.Fst ~ altitude,
      data = surface,
      conductance_model = factory,
      measurement_model = leastsquares,
      theta = c(altitude = 0, "sigma.altitude" = 5),
      optimizer = "bfgs",
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
    ),
    "must lie within the conductance-model bounds"
  )

  expect_error(
    terradish(
      melip.Fst ~ altitude,
      data = surface,
      conductance_model = factory,
      measurement_model = leastsquares,
      approximation = "coarse_raster",
      optimizer = "bfgs",
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
    ),
    "not currently supported"
  )
})
