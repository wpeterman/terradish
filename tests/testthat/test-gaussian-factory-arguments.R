# Regression test for gaussian_smoothed_loglinear_conductance() argument
# binding. The factory and its coarse-raster rebuild helper read
# scale_vars, standardize, sigma_lower, sigma_upper and
# sigma_conversion_factor after the constructor returns. Before 0.0.60 those
# arguments stayed unevaluated promises, so factories built in a for loop all
# read the loop variable's final value: every factory built with
# scale_vars = v smoothed the last v. These checks only instantiate the
# factories (no fitting), so they run in well under a second.

test_that("Gaussian factories built in a for loop keep their own arguments", {
  dat <- melip_fixture(1:6)
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8,
                                 saveStack = TRUE)
  rhs <- ~ altitude + forestcover
  vars <- c("altitude", "forestcover")
  cell <- mean(terra::res(dat$covariates))

  # Each check loop uses a different variable from its build loop. Reusing the
  # build variable would re-bind it to the expected value just before each
  # factory is first called and hide the bug.

  # scale_vars: each factory smooths the layer it was built with
  by_var <- list()
  for (v in vars)
    by_var[[v]] <- gaussian_smoothed_loglinear_conductance(surface, scale_vars = v)
  for (want in vars) {
    model <- by_var[[want]](rhs, surface$x)
    expect_identical(attr(model, "scale_vars"), want)
    expect_identical(grep("^sigma\\.", names(attr(model, "default")), value = TRUE),
                     paste0("sigma.", want))
  }

  # sigma_upper: each factory keeps its own bound (map units)
  uppers <- c(2, 3) * cell
  by_upper <- list()
  for (u in uppers)
    by_upper[[format(u)]] <- gaussian_smoothed_loglinear_conductance(
      surface, scale_vars = "altitude", sigma_lower = 0.5 * cell, sigma_upper = u)
  for (want in uppers) {
    info <- attr(by_upper[[format(want)]](rhs, surface$x), "gaussian_scale_info")
    expect_equal(unname(info$upper[["altitude"]]), want)
  }

  # sigma_lower: each factory keeps its own bound (map units)
  lowers <- c(0.5, 1) * cell
  by_lower <- list()
  for (l in lowers)
    by_lower[[format(l)]] <- gaussian_smoothed_loglinear_conductance(
      surface, scale_vars = "altitude", sigma_lower = l, sigma_upper = 3 * cell)
  for (want in lowers) {
    info <- attr(by_lower[[format(want)]](rhs, surface$x), "gaussian_scale_info")
    expect_equal(unname(info$lower[["altitude"]]), want)
  }

  # standardize: each factory keeps its own setting
  by_std <- list()
  for (s in c(TRUE, FALSE))
    by_std[[as.character(s)]] <- gaussian_smoothed_loglinear_conductance(
      surface, scale_vars = "altitude", standardize = s)
  for (want in c(TRUE, FALSE))
    expect_identical(attr(by_std[[as.character(want)]](rhs, surface$x),
                          "gaussian_scale_info")$standardize, want)

  # sigma_conversion_factor: each factory keeps its own conversion
  factors <- c(1, 2)
  by_factor <- list()
  for (k in factors)
    by_factor[[format(k)]] <- gaussian_smoothed_loglinear_conductance(
      surface, scale_vars = "altitude", sigma_conversion_factor = k)
  for (want in factors) {
    info <- attr(by_factor[[format(want)]](rhs, surface$x), "gaussian_scale_info")
    expect_equal(unname(info$conversion[["altitude"]]), want)
  }

  # the coarse-raster rebuild helper reads the same arguments, through a fresh
  # set of factories that are never called directly
  fresh <- list()
  for (v in vars)
    fresh[[v]] <- gaussian_smoothed_loglinear_conductance(surface, scale_vars = v)
  for (want in vars) {
    rebuilt <- attr(fresh[[want]], "rebuild_for_surface")(rhs, surface)
    expect_identical(attr(rebuilt, "scale_vars"), want)
  }
})
