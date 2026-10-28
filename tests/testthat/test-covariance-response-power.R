# Simulation-based power regressions run in CI with NOT_CRAN=true.
skip_on_cran()

test_that("power summaries exclude failed outer and nuisance fits from selection", {
  good <- structure(list(mle = list(theta = c(x = .2)), loglik = -10,
    aic = 24, df = 2, convergence = list(code = 0L),
    fit = list(subproblem = list(convergence = 0L)),
    submodels = list(f = function(theta) list(conductance = 1:3))), class = "terradish")
  row <- function(fit, label) .covariance_power_result_row(
    "test", "spacefill", 6, 1, 1, label, "OK", 0, 6, 15, fit, 1:3, .8)
  for (part in c("outer", "nuisance", "missing")) {
    bad <- good
    # A higher apparent likelihood must not make an incomplete fit selectable.
    bad$loglik <- 100; bad$aic <- -196
    if (part == "outer") bad$convergence$code <- 1L
    if (part == "nuisance") bad$fit$subproblem$convergence <- 2L
    if (part == "missing") bad$fit$subproblem <- NULL
    result <- row(bad, "bad")
    expect_false(result$fit_ok)
    expect_match(result$status, "^NONCONVERGED:")
    expect_true(all(is.na(result[c("loglik", "AIC", "AICc")])))
    expect_false(result$conductance_recovered)
    ranked <- .covariance_power_mark_selection(rbind(row(good, "good"), result))
    expect_identical(ranked$selected_AIC, c(TRUE, FALSE))
    expect_identical(ranked$selected_AICc, c(TRUE, FALSE))
  }
  expect_true(row(good, "good")$fit_ok)
  expect_true(row(good, "good")$conductance_recovered)
})

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
    sigma_upper = 1,
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
