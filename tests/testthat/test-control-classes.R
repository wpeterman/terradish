# NewtonRaphsonControl() and HagerZhangControl() now return classed lists, as
# ArmijoControl() always has. The class is a label only: the elements, their
# order, and every fitted number must be unchanged.

test_that("control constructors return classed lists with unchanged elements", {
  newton <- NewtonRaphsonControl()
  expect_s3_class(newton, "terradish_newton_raphson_control")
  expect_true(is.list(newton))
  expect_named(newton, c("maxit", "ctol", "etol", "ftol", "verbose", "eps",
                         "del", "ls.control"))

  hager_zhang <- HagerZhangControl()
  expect_s3_class(hager_zhang, "terradish_hager_zhang_control")
  expect_true(is.list(hager_zhang))
  expect_named(hager_zhang, c("delta", "sigma", "alphamax", "rho", "epsilon",
                              "gamma", "linesearchmax", "psi3", "c",
                              "verbose"))

  # The default line search is the classed Hager-Zhang control, and an Armijo
  # control is still recognized by its own class.
  expect_s3_class(newton$ls.control, "terradish_hager_zhang_control")
  armijo <- NewtonRaphsonControl(ls.control = ArmijoControl())
  expect_s3_class(armijo$ls.control, "terradish_armijo_control")
  expect_false(inherits(newton$ls.control, "terradish_armijo_control"))

  # Arguments are stored unchanged.
  custom <- NewtonRaphsonControl(maxit = 7, ctol = 1e-4,
                                 ls.control = HagerZhangControl(sigma = 0.5))
  expect_identical(custom$maxit, 7)
  expect_identical(custom$ctol, 1e-4)
  expect_identical(custom$ls.control$sigma, 0.5)
})

test_that("classed and unclassed controls give identical fits", {
  dat <- melip_fixture(1:8, aggregate_factor = 4)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords)

  classed <- NewtonRaphsonControl(maxit = 50, verbose = FALSE)
  # A plain list, as the constructors returned before, must still be accepted.
  plain <- unclass(classed)
  plain$ls.control <- unclass(plain$ls.control)
  expect_identical(class(plain), "list")

  for (optimizer in c("newton", "bfgs")) {
    fit_classed <- terradish(melip.Fst ~ altitude + forestcover, surface,
                             measurement_model = mlpe, optimizer = optimizer,
                             control = classed)
    fit_plain <- terradish(melip.Fst ~ altitude + forestcover, surface,
                           measurement_model = mlpe, optimizer = optimizer,
                           control = plain)
    expect_identical(coef(fit_classed), coef(fit_plain))
    expect_identical(fit_classed$loglik, fit_plain$loglik)
    expect_identical(vcov(fit_classed), vcov(fit_plain))
    expect_identical(fit_classed$fit$phi, fit_plain$fit$phi)
  }
})

test_that("settings assessment still labels the line search from the class", {
  expect_identical(
    .terradish_assessment_line_search_name(HagerZhangControl()),
    "hager_zhang")
  expect_identical(
    .terradish_assessment_line_search_name(ArmijoControl()),
    "armijo")
})
