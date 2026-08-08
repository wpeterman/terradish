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

# A synthetic lattice with a known conductance truth. Small enough to fit
# repeatedly, and well enough determined that a large tau2 makes the nuisance
# subproblem singular, which is the failure this guards against.
hierarchical_fixture <- function() {
  r  <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12,
                    ymin = 0, ymax = 12)
  gx <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  gy <- terra::yFromCell(r, seq_len(terra::ncell(r)))
  covariates <- c(terra::setValues(r, scale(gx)[, 1]),
                  terra::setValues(r, scale(gy)[, 1]))
  names(covariates) <- c("altitude", "forestcover")
  coords <- terra::xyFromCell(r, c(1, 12, 66, 79, 133, 144, 40, 105))
  surface <- conductance_surface(covariates, coords, directions = 8,
                                 saveStack = TRUE)
  E <- terradish_algorithm(
    loglinear_conductance(~ altitude + forestcover, surface$x),
    leastsquares, surface, S = diag(nrow(coords)), theta = c(0.5, -0.4),
    objective = FALSE, gradient = FALSE, hessian = FALSE,
    partial = FALSE)$covariance
  E <- as.matrix(E)
  S <- outer(diag(E), diag(E), "+") - 2 * E
  diag(S) <- 0
  list(surface = surface, S = S, truth = c(altitude = 0.5, forestcover = -0.4))
}

test_that("a failing tau2 grid point is skipped, recorded, and reported", {
  # tau2 = 10 makes the nuisance Hessian singular from the warm start carried
  # in from tau2 = 1. Before the fix this aborted the whole fit. One verbose
  # fit exercises every part of the contract, so the fixture is built once.
  fx <- hierarchical_fixture()
  S <- fx$S
  warnings_seen <- character()
  msgs <- character()

  fit <- withCallingHandlers(
    terradish_hierarchical(S ~ altitude + forestcover, data = fx$surface,
                           measurement_model = leastsquares,
                           field_resolution = 3L, tau2_grid = c(0.1, 1, 10),
                           verbose = TRUE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })

  # the fit survives
  expect_s3_class(fit, "terradish_hierarchical")
  expect_true(is.finite(fit$tau2))

  # the failed point is recorded as NA rather than dropped from the table
  expect_equal(nrow(fit$tau2_selection), 3L)
  expect_true(any(is.na(fit$tau2_selection$logML)))
  expect_true(all(is.finite(fit$tau2_selection$logML[1:2])))

  # tau2 is chosen from the points that worked
  usable <- fit$tau2_selection$tau2[is.finite(fit$tau2_selection$logML)]
  expect_true(fit$tau2 %in% usable)

  # the covariate effects are still recovered
  expect_equal(unname(coef(fit)), unname(fx$truth), tolerance = 0.05)

  # the skip is surfaced, not silent, and the grid-edge warning describes the
  # usable range rather than the requested one
  expect_true(any(grepl("were skipped because the fit failed", warnings_seen)))
  expect_true(any(grepl("edge of the usable", warnings_seen)))
  expect_true(any(grepl("skipped", msgs)))
  expect_true(any(grepl("logML", msgs)))
})

test_that("terradish_hierarchical is quiet by default", {
  fx <- hierarchical_fixture()
  S <- fx$S
  expect_no_message(
    suppressWarnings(
      terradish_hierarchical(S ~ altitude + forestcover, data = fx$surface,
                             measurement_model = leastsquares,
                             field_resolution = 3L, tau2 = 0.1)))
})
