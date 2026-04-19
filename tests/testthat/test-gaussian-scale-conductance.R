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
  info <- attr(conductance_model, "gaussian_scale_info")
  expect_equal(info$lower[["altitude"]], 0.5)
  expect_equal(info$upper[["altitude"]], 2)
  expect_equal(
    attr(conductance_model, "lower")[["sigma.altitude"]],
    0.5 / info$conversion[["altitude"]]
  )
  expect_equal(
    attr(conductance_model, "upper")[["sigma.altitude"]],
    2 / info$conversion[["altitude"]]
  )

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
  expect_true(all(c("covariate", "sigma", "sigma_internal",
                    "sigma_conversion", "sigma_cells_x", "radial_95",
                    "sigma_km", "radial_95_km") %in% names(sigma_table)))
  expect_equal(sigma_table$covariate, "altitude")
  expect_equal(
    unname(coef(fit)[["sigma.altitude"]]),
    unname(fit$mle$theta_internal[["sigma.altitude"]]) *
      sigma_table$sigma_conversion[1],
    tolerance = 1e-10
  )
})

test_that("gaussian sigma conversion factor keeps map-unit reporting", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["altitude"]]
  names(covariates) <- "altitude"
  surface <- conductance_surface(covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)

  fit <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude,
      data = surface,
      conductance_model = gaussian_smoothed_loglinear_conductance(
        surface,
        sigma_lower = 0.5,
        sigma_upper = 2,
        sigma_conversion_factor = 0.25
      ),
      measurement_model = leastsquares,
      optimizer = "bfgs",
      leverage = FALSE,
      control = NewtonRaphsonControl(maxit = 6, verbose = FALSE)
    )
  )

  expect_equal(
    unname(coef(fit)[["sigma.altitude"]]),
    unname(fit$mle$theta_internal[["sigma.altitude"]]) * 0.25,
    tolerance = 1e-10
  )

  sigma_table <- suppressWarnings(gaussian_scale_summary(fit))
  expect_equal(sigma_table$sigma_conversion, 0.25)
  expect_equal(sigma_table$sigma_internal, sigma_table$sigma / 0.25)
})

test_that("gaussian marginal plots work on the fitted smoothed covariate scale", {
  dat <- melip_fixture(1:6)
  melip.Fst <- dat$melip.Fst
  covariates <- dat$covariates[["forestcover"]]
  names(covariates) <- "forestcover"
  surface <- conductance_surface(covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)

  fit <- suppressWarnings(
    terradish(
      melip.Fst ~ forestcover,
      data = surface,
      conductance_model = gaussian_smoothed_loglinear_conductance(
        surface,
        sigma_lower = 0.5,
        sigma_upper = 2
      ),
      measurement_model = leastsquares,
      optimizer = "bfgs",
      leverage = FALSE,
      control = NewtonRaphsonControl(maxit = 6, verbose = FALSE)
    )
  )

  marg <- plot(fit, type = "marginal", data = surface, n = 12)
  expect_s3_class(marg, "ggplot")
  expect_true(all(c("covariate", "x", "est", "lower", "upper") %in%
                    names(marg$data)))
  expect_true(any(grepl("smoothed", levels(marg$data$covariate), fixed = TRUE)))
  expect_true(all(marg$data$lower <= marg$data$est))
  expect_true(all(marg$data$est <= marg$data$upper))

  marg_response <- plot(fit, type = "marginal_response", data = surface, n = 8)
  expect_s3_class(marg_response, "ggplot")
  expect_equal(marg_response$scales$get_scales("y")$name,
               "Predicted genetic distance")
  expect_true(all(is.finite(marg_response$data$est)))

  sigma_plot <- plot(fit, type = "sigma", n = 40)
  expect_s3_class(sigma_plot, "ggplot")
  expect_true(all(c("covariate", "distance", "weight") %in% names(sigma_plot$data)))
  expect_true(all(sigma_plot$data$distance >= 0))
  expect_true(all(sigma_plot$data$weight >= 0))

  sigma_plot_km <- plot(
    fit,
    type = "sigma",
    n = 40,
    distance_per_map_unit = 0.001,
    distance_unit = "km"
  )
  expect_s3_class(sigma_plot_km, "ggplot")
  expect_equal(sigma_plot_km$labels$x, "Distance (km)")
})

test_that("gaussian scale-aware conductance rejects factor-valued raster layers", {
  dat <- melip_fixture(1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)
  surface$x$forestcover <- factor(ifelse(surface$x$forestcover > 0, "high", "low"))

  expect_error(
    gaussian_smoothed_loglinear_conductance(surface)(
      ~ altitude + forestcover,
      surface$x
    ),
    "Factor-valued raster layers are not supported"
  )
})

test_that("gaussian scale-aware conductance validates sigma starts and supports coarse raster starts", {
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

  fit_coarse <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude,
      data = surface,
      conductance_model = factory,
      measurement_model = leastsquares,
      approximation = "coarse_raster",
      approximation_control = list(factor = 2L, exact_refine = TRUE),
      optimizer = "bfgs",
      leverage = FALSE,
      control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
    )
  )

  expect_s3_class(fit_coarse, "terradish")
  expect_true(isTRUE(fit_coarse$approximation$used))
  expect_equal(fit_coarse$approximation$type, "coarse_raster")
  expect_true(all(c("altitude", "sigma.altitude") %in% names(coef(fit_coarse))))
})
