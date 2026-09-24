# Two behavioral contracts that are easy to regress silently:
#   1. fitting is quiet by default, and the trace is a suppressible message()
#      rather than cat() output;
#   2. a tau2 grid point whose fit fails is skipped, not fatal.

test_that("terradish() is quiet by default and verbose = TRUE turns the trace on", {
  res <- fit_fixture(keep = 1:8)
  surface <- res$surface
  melip.Fst <- res$data$melip.Fst

  fit_quietly <- function(...)
    suppressWarnings(
      terradish(melip.Fst ~ altitude + forestcover, data = surface,
                conductance_model = loglinear_conductance,
                measurement_model = leastsquares,
                control = NewtonRaphsonControl(maxit = 3), ...))

  # nothing on stdout or stderr with the defaults
  expect_silent(fit_default <- fit_quietly())

  # verbose = TRUE emits messages, and messages specifically (not cat output),
  # so users can silence them
  expect_message(fit_loud <- fit_quietly(verbose = TRUE))
  expect_silent(suppressMessages(fit_quietly(verbose = TRUE)))

  # the switch must not change the answer
  expect_equal(coef(fit_default), coef(fit_loud))
  expect_equal(fit_default$loglik, fit_loud$loglik)
})

test_that("an explicit verbose overrides control, and control alone still works", {
  res <- fit_fixture(keep = 1:8)
  surface <- res$surface
  melip.Fst <- res$data$melip.Fst

  loud_control <- NewtonRaphsonControl(maxit = 3, verbose = TRUE)

  # control on its own is honored
  expect_message(
    suppressWarnings(
      terradish(melip.Fst ~ altitude + forestcover, data = surface,
                conductance_model = loglinear_conductance,
                measurement_model = leastsquares, control = loud_control)))

  # verbose = FALSE wins over a loud control object
  expect_silent(
    suppressWarnings(
      terradish(melip.Fst ~ altitude + forestcover, data = surface,
                conductance_model = loglinear_conductance,
                measurement_model = leastsquares, control = loud_control,
                verbose = FALSE)))

  expect_error(
    terradish(melip.Fst ~ altitude + forestcover, data = surface,
              conductance_model = loglinear_conductance,
              measurement_model = leastsquares, verbose = "yes"),
    "must be TRUE or FALSE")
})

test_that("optimizer control objects default to quiet", {
  expect_false(NewtonRaphsonControl()$verbose)
  expect_false(HagerZhangControl()$verbose)
  expect_false(ArmijoControl()$verbose)
})
