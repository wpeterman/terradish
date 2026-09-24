prediction_fixture <- function(factory = loglinear_conductance, formula = ~x) {
  g <- robustness_surface()
  if (identical(factory, "gaussian")) factory <-
    gaussian_smoothed_loglinear_conductance(g, sigma_lower = 0.3, sigma_upper = 1)
  f <- factory(formula, g$x)
  theta <- attr(f, "default")
  theta[] <- 0.2
  if (any(grepl("sigma", names(theta)))) theta[grepl("sigma", names(theta))] <- 0.7
  fit <- structure(list(formula = formula, fit = list(boundary = FALSE),
    submodels = list(f_internal = f, f_factory = factory, f = f, g = leastsquares),
    mle = list(theta = theta, hessian = -diag(length(theta)))), class = "terradish")
  list(g = g, f = f, factory = factory, fit = fit, theta = theta)
}

test_that("new log-linear landscapes retain terms and fitted transformations", {
  for (formula in list(~x, ~I(x^2), ~poly(x, 2), ~scale(x))) {
    z <- prediction_fixture(formula = formula)
    same <- conductance(z$g, z$fit)
    expect_equal(c(terra::values(same[[1]])), z$f(z$theta)$conductance)
    shifted <- z$g$stack + 0.4
    g2 <- conductance_surface(shifted, z$g$vertex_coordinates[z$g$demes, ])
    predicted <- conductance(g2, z$fit)
    if (identical(formula, ~x)) expected <- exp(g2$x$x * z$theta)
    else {
      spec <- attr(z$f, "prediction_spec")
      frame <- model.frame(spec$terms, g2$x)
      X <- model.matrix(spec$terms, frame)[, names(z$theta), drop = FALSE]
      expected <- exp(c(X %*% z$theta))
    }
    expect_equal(c(terra::values(predicted[[1]])), c(expected), tolerance = 1e-11)
  }
})

test_that("new spline landscapes retain training knots and centering", {
  z <- prediction_fixture(smooth_loglinear_conductance, ~s(x, df = 3))
  expect_equal(c(terra::values(conductance(z$g, z$fit)[[1]])), z$f(z$theta)$conductance)
  shifted <- z$g$stack + 0.1
  g2 <- conductance_surface(shifted, z$g$vertex_coordinates[z$g$demes, ])
  spec <- attr(z$f, "smooth_loglinear_info")$smooth_specs[[1]]
  B <- splines::ns(g2$x$x, knots = spec$knots, Boundary.knots = spec$Boundary.knots)
  expected <- exp(c(sweep(B, 2, spec$centers, "-") %*% z$theta))
  expect_equal(c(terra::values(conductance(g2, z$fit)[[1]])), expected, tolerance = 1e-11)
})

test_that("mixed spline models retain transformations of parametric terms", {
  x <- data.frame(x = seq(-2, 2, length.out = 60), y = cos(1:60))
  formula <- ~scale(y) + s(x, df = 3)
  f <- smooth_loglinear_conductance(formula, x)
  theta <- setNames(rep(0.2, 4), names(attr(f, "default")))
  new <- transform(x, y = y + 0.5)
  rebuilt <- attr(f, "plot_factory")(formula, new)
  expect_equal(log(rebuilt(theta)$conductance) - log(f(theta)$conductance),
    rep(0.2 * 0.5 / sd(x$y), nrow(x)), tolerance = 1e-11)
})

test_that("Gaussian prediction uses training post-smoothing scaling and derivatives", {
  z <- prediction_fixture("gaussian")
  expect_equal(c(terra::values(conductance(z$g, z$fit)[[1]])), z$f(z$theta)$conductance)
  shifted <- z$g$stack + 0.6
  g2 <- conductance_surface(shifted, z$g$vertex_coordinates[z$g$demes, ])
  cells <- terra::cellFromXY(z$g$stack, z$g$vertex_coordinates)
  rc <- terra::rowColFromCell(z$g$stack, cells)
  train <- .gaussian_scale_layer_values(.gaussian_scale_prepare_layer(z$g$stack, rc), 0.7, FALSE)$value
  test <- .gaussian_scale_layer_values(.gaussian_scale_prepare_layer(shifted, rc), 0.7, FALSE)$value
  expected <- exp(0.2 * (test - mean(train)) / sd(train))
  expect_equal(c(terra::values(conductance(g2, z$fit)[[1]])), expected, tolerance = 1e-10)
  rebuilt <- attr(z$factory, "predict_for_surface")(~x, g2, z$f)
  numerical <- numDeriv::jacobian(function(theta) rebuilt(theta)$conductance, z$theta)
  analytic <- vapply(1:2, function(i) rebuilt(z$theta)$df__dtheta(i), numeric(nrow(g2$x)))
  expect_equal(unname(analytic), unname(numerical), tolerance = 1e-7)
  expect_equal(z$f(z$theta)$conductance, prediction_fixture("gaussian")$f(z$theta)$conductance)
})
