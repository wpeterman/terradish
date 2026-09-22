test_that("slim_terradish reduces storage and retains inference methods", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))
  fit <- fx$fit
  slim <- slim_terradish(fit)

  expect_s3_class(slim, "terradish")
  expect_null(slim$submodels)
  expect_null(slim$leverage)
  expect_lt(as.numeric(object.size(slim)), as.numeric(object.size(fit)))
  expect_equal(coef(slim), coef(fit))
  expect_equal(as.numeric(logLik(slim)), as.numeric(logLik(fit)))
  expect_equal(fitted(slim), fitted(fit))
  expect_s3_class(summary(slim), "summary.terradish")
  expect_error(conductance(fx$surface, slim), "slim_terradish")
  expect_error(plot(slim, type = "surface", data = fx$surface),
               "slim_terradish")
})

test_that("terradish slim option and response dropping are explicit", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 1, verbose = FALSE))
  melip.Fst <- fx$data$melip.Fst
  fit <- suppressWarnings(terradish(
    melip.Fst ~ altitude, fx$surface,
    measurement_model = leastsquares,
    control = NewtonRaphsonControl(maxit = 1, verbose = FALSE),
    slim = TRUE
  ))
  expect_true(isTRUE(fit$storage$slim))
  expect_null(fit$submodels)

  no_response <- slim_terradish(fx$fit, drop = "response")
  expect_null(no_response$fit$response)
  expect_null(no_response$comparison$response)
  expect_error(slim_terradish(fx$fit, drop = "unknown"), "Unknown component")
})
