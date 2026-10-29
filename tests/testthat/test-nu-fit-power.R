test_that("overstated fit information reduces coverage without changing the truth", {
  # Thirty simulation replicates are retained in the extended suite.
  skip_on_cran()
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
  # nu_fit = 20 * nu scales the log-likelihood by 20, so both fits share one
  # maximum. A replicate with two likelihood peaks can still send the local
  # optimizer to different peaks across platforms, so compare only replicates
  # whose fits reached the same scaled log-likelihood.
  loglik_at <- function(x)
    x$results$loglik[match(x$parameter_results$simulation, x$results$simulation)]
  same_peak <- keep &
    abs(loglik_at(overstated) - 20 * loglik_at(matched)) < 1e-3
  expect_gt(sum(same_peak), 20)
  expect_equal(matched$parameter_results$estimate[same_peak],
    overstated$parameter_results$estimate[same_peak], tolerance = 1e-5)
  expect_equal(matched$parameter_results$se[same_peak] / sqrt(20),
    overstated$parameter_results$se[same_peak], tolerance = 1e-5)
})
