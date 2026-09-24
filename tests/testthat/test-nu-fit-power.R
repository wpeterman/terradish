test_that("overstated fit information reduces coverage without changing the truth", {
  g <- robustness_surface()
  args <- list(theta = c(x = 0.6), formula = ~x, data = g,
    sample_sizes = 6, strategies = "spacefill", nu = 50, nsim = 30,
    tau = 0.7, sigma = 0.08, seed = 905,
    control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))
  matched <- do.call(covariance_response_power, args)
  overstated <- do.call(covariance_response_power, c(args, list(nu_fit = 1000)))
  expect_equal(matched$settings$nu_fit, 50)
  expect_equal(overstated$settings$nu, 50)
  expect_equal(overstated$settings$nu_fit, 1000)
  expect_lt(overstated$parameter_summary$coverage, matched$parameter_summary$coverage)
  keep <- matched$parameter_results$fit_ok & overstated$parameter_results$fit_ok
  expect_gt(sum(keep), 20)
  expect_equal(matched$parameter_results$estimate[keep],
    overstated$parameter_results$estimate[keep], tolerance = 1e-5)
  expect_equal(matched$parameter_results$se[keep] / sqrt(20),
    overstated$parameter_results$se[keep], tolerance = 1e-5)
})
