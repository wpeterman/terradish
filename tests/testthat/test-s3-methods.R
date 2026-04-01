test_that("S3 methods for radish objects return consistent outputs", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  fit <- fx$fit

  expect_gt(length(capture.output(print(fit))), 0)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.radish")
  expect_gt(length(capture.output(print(sm))), 0)

  expect_equal(coef(fit), fit$mle$theta)
  expect_true(is.matrix(fitted(fit, type = "response")))
  expect_true(is.matrix(fitted(fit, type = "distance")))
  expect_true(is.matrix(fitted(fit, type = "covariance")))
  expect_true(is.matrix(resid(fit)))
  expect_equal(dim(simulate(fit, nsim = 2)), c(nrow(fit$fit$response), ncol(fit$fit$response), 2))
  expect_s3_class(logLik(fit), "logLik")
  expect_equal(AIC(fit), fit$aic)
})

test_that("anova compares fitted terradish models", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  surface <- fx$surface
  melip.Fst <- fx$data$melip.Fst

  fit1 <- suppressWarnings(
    radish(melip.Fst ~ altitude,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )
  fit2 <- suppressWarnings(
    radish(melip.Fst ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = leastsquares,
           control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  )

  tab <- anova(fit1, fit2)
  expect_true(is.matrix(tab) || is.data.frame(tab))
})
