# Two options were offered but never implemented. simulate() no longer lists
# "parametric", and terradish(validate = TRUE) says plainly that the argument
# is deprecated and ignored.

test_that("simulate() offers only the implemented permutation method", {
  fx <- fit_fixture(keep = 1:8, aggregate_factor = 4)

  expect_identical(eval(formals(simulate.terradish)$method), "permutation")

  sim <- simulate(fx$fit, seed = 1)
  expect_equal(dim(sim), dim(fitted(fx$fit)))
  expect_identical(sim, simulate(fx$fit, seed = 1, method = "permutation"))

  expect_error(simulate(fx$fit, method = "parametric"),
               "`method` must be \"permutation\"")
  expect_error(simulate(fx$fit, method = c("permutation", "parametric")),
               "`method` must be \"permutation\"")
})

test_that("terradish(validate = TRUE) warns that the argument is ignored", {
  dat <- melip_fixture(1:8, aggregate_factor = 4)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords)
  control <- NewtonRaphsonControl(maxit = 50, verbose = FALSE)

  expect_warning(
    validated <- terradish(melip.Fst ~ altitude + forestcover, surface,
                           measurement_model = leastsquares,
                           control = control, validate = TRUE),
    "`validate` argument of terradish\\(\\) is deprecated and ignored"
  )
  plain <- terradish(melip.Fst ~ altitude + forestcover, surface,
                     measurement_model = leastsquares, control = control)

  # The argument changes nothing about the fit.
  expect_identical(coef(validated), coef(plain))
  expect_identical(validated$loglik, plain$loglik)
  expect_identical(vcov(validated), vcov(plain))
  expect_null(validated$leverage$validate)
})
