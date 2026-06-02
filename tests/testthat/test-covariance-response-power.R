test_that("covariance_response_power summarizes spline recovery scenarios", {
  dat <- melip_fixture(keep = 1:7)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  true_formula <- ~ forestcover + s(altitude, df = 3)
  true_model <- smooth_loglinear_conductance(true_formula, surface$x)
  theta_true <- attr(true_model, "default")
  theta_true["forestcover"] <- -0.1
  theta_true[grepl("^s\\(altitude\\)", names(theta_true))] <-
    c(-0.35, 0.5, -0.25)

  power <- covariance_response_power(
    theta = theta_true,
    formula = true_formula,
    data = surface,
    sample_sizes = 6,
    strategies = "spacefill",
    conductance_model = smooth_loglinear_conductance,
    fit_models = list(
      spline = list(formula = true_formula,
                    conductance_model = smooth_loglinear_conductance)
    ),
    tau = 0.7,
    sigma = 0.12,
    nu = 50,
    nsim = 1,
    seed = 101,
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
  )

  expect_s3_class(power, "terradish_covariance_power")
  expect_equal(nrow(power$results), 1)
  expect_equal(power$results$model, "spline")
  expect_equal(power$results$sample_size, 6)
  expect_true(all(c("conductance_power", "all_parameter_power") %in%
                    names(power$summary)))
  expect_true(nrow(power$parameter_summary) >= 1)
  expect_true(all(power$parameter_summary$model == "spline"))
})

test_that("covariance_response_power accepts Gaussian theta on external scale", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)
  gaussian_factory <- gaussian_smoothed_loglinear_conductance(
    surface,
    scale_vars = "altitude",
    sigma_lower = 0.5,
    sigma_upper = 4,
    sigma_conversion_factor = 2
  )
  gaussian_model <- gaussian_factory(~ altitude, surface$x)
  theta_external <- attr(gaussian_model, "default")
  theta_external["altitude"] <- 0.15
  theta_external["sigma.altitude"] <- theta_external["sigma.altitude"] * 2

  power <- covariance_response_power(
    theta = theta_external,
    formula = ~ altitude,
    data = surface,
    sample_sizes = 5,
    strategies = "spacefill",
    conductance_model = gaussian_factory,
    tau = 0.7,
    sigma = 0.12,
    nu = 40,
    nsim = 1,
    seed = 102,
    optimizer = "bfgs",
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
  )

  expect_equal(power$settings$theta["sigma.altitude"],
               power$settings$theta_internal["sigma.altitude"] * 2)
  expect_equal(power$results$model, "truth")
  expect_equal(power$summary$sample_size, 5)
})

test_that("covariance_response_power counts boundary fits as non-detections", {
  dat <- melip_fixture(keep = 1:6)
  surface <- conductance_surface(dat$covariates[["altitude"]], dat$coords,
                                 directions = 8)

  power <- covariance_response_power(
    theta = c(altitude = 0.2),
    formula = ~ altitude,
    data = surface,
    sample_sizes = 6,
    strategies = "spacefill",
    tau = 0,
    sigma = 1,
    nu = 100000,
    nsim = 1,
    seed = 103,
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE)
  )

  expect_false(power$results$fit_ok)
  expect_equal(nrow(power$parameter_summary), 1)
  expect_equal(power$parameter_summary$fit_rate, 0)
  expect_equal(power$parameter_summary$power, 0)
})
