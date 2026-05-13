test_that("smooth_loglinear_conductance expands smooth terms", {
  x <- data.frame(
    altitude = seq(-1, 1, length.out = 12),
    forestcover = seq(0.1, 0.9, length.out = 12)
  )

  model <- smooth_loglinear_conductance(
    ~ forestcover + s(altitude, df = 3),
    x
  )
  theta <- attr(model, "default")

  expect_s3_class(model, "terradish_conductance_model")
  expect_length(theta, 4)
  expect_named(theta, c("forestcover", paste0("s(altitude).", 1:3)))

  fit <- model(c(forestcover = 0.2, `s(altitude).1` = -0.1,
                 `s(altitude).2` = 0.3, `s(altitude).3` = -0.2))

  expect_length(fit$conductance, nrow(x))
  expect_true(all(fit$conductance > 0))
  expect_equal(dim(fit$df__dtheta_matrix), c(nrow(x), length(theta)))
})

test_that("smooth_loglinear_conductance derivatives match finite differences", {
  x <- data.frame(
    altitude = seq(-1, 1, length.out = 12),
    forestcover = seq(0.1, 0.9, length.out = 12)
  )
  model <- smooth_loglinear_conductance(
    ~ forestcover + s(altitude, df = 3),
    x
  )
  theta <- c(forestcover = 0.2, `s(altitude).1` = -0.1,
             `s(altitude).2` = 0.3, `s(altitude).3` = -0.2)
  fit <- model(theta)

  step <- 1e-6
  for (k in seq_along(theta)) {
    delta <- rep(0, length(theta))
    delta[k] <- step
    finite_difference <- (
      model(theta + delta)$conductance -
        model(theta - delta)$conductance
    ) / (2 * step)

    expect_equal(fit$df__dtheta(k), finite_difference, tolerance = 1e-6)
  }
})

test_that("smooth_loglinear_conductance supports k alias and B-splines", {
  x <- data.frame(altitude = seq(0, 10, length.out = 15))

  model <- smooth_loglinear_conductance(
    ~ s(altitude, k = 4, basis = "bs", degree = 2),
    x
  )
  theta <- attr(model, "default")

  expect_length(theta, 4)
  expect_named(theta, paste0("s(altitude).", 1:4))
  expect_true(all(is.finite(model(theta)$conductance)))
})

test_that("smooth_loglinear_conductance rejects non-finite conductance values", {
  x <- data.frame(
    altitude = seq(0, 10, length.out = 15),
    forestcover = seq(10, 150, length.out = 15)
  )
  model <- smooth_loglinear_conductance(
    ~ forestcover + s(altitude, df = 3),
    x
  )
  theta <- attr(model, "default")
  theta["forestcover"] <- 1e308

  expect_error(
    model(theta),
    "non-finite conductance values"
  )
})

test_that("smooth_loglinear_conductance carries plotting metadata", {
  x <- data.frame(altitude = seq(0, 10, length.out = 15))
  model <- smooth_loglinear_conductance(
    ~ s(altitude, k = 4, basis = "bs", degree = 2),
    x
  )
  factory <- attr(model, "plot_factory")
  rebuilt <- factory(~ s(altitude), x)

  expect_s3_class(factory, "terradish_conductance_model_factory")
  expect_identical(attr(factory, "link"), "log")
  expect_identical(attr(model, "link"), "log")
  expect_equal(names(attr(rebuilt, "default")), paste0("s(altitude).", 1:4))
})

test_that("plot methods resolve smooth conductance models automatically", {
  fx <- fit_fixture(keep = 1:8,
                    control = NewtonRaphsonControl(maxit = 2,
                                                   verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst

  fit <- suppressWarnings(
    terradish(melip.Fst ~ s(altitude, df = 3),
              data = surface,
              conductance_model = smooth_loglinear_conductance,
              measurement_model = leastsquares,
              control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  expect_true(grepl("^s\\(altitude\\)\\.", names(coef(fit))[1]))

  marg <- plot(fit, type = "marginal", data = surface, n = 10)
  expect_true(inherits(marg, "ggplot") || is.list(marg))
  marg_data <- if (inherits(marg, "ggplot")) marg$data else do.call(rbind, lapply(marg, function(p) p$data))
  expect_true("altitude (original scale)" %in%
                unique(as.character(marg_data$covariate)))
  expect_true(all(c("x", "est", "lower", "upper") %in% names(marg_data)))
  expect_true(all(is.finite(marg_data$est)))

  marg_response <- plot(fit, type = "marginal_response", data = surface, n = 5)
  expect_true(inherits(marg_response, "ggplot") || is.list(marg_response))
  marg_response_data <- if (inherits(marg_response, "ggplot")) marg_response$data else do.call(rbind, lapply(marg_response, function(p) p$data))
  expect_true("altitude (original scale)" %in%
                unique(as.character(marg_response_data$covariate)))
  expect_true(all(is.finite(marg_response_data$est)))
})

test_that("smooth_loglinear_conductance rejects unsupported smooth terms", {
  x <- data.frame(altitude = seq(0, 10, length.out = 15))

  expect_error(
    smooth_loglinear_conductance(~ s(I(altitude^2), df = 3), x),
    "one raw column name"
  )
  expect_error(
    smooth_loglinear_conductance(~ s(missing_layer, df = 3), x),
    "not found"
  )
})
