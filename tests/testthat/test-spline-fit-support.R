test_that("focal spline support and exact coarse refinement reuse the fitted basis", {
  data(melip, package = "terradish", envir = environment())
  rasters <- c(terra::unwrap(melip.altitude), terra::unwrap(melip.forestcover))
  names(rasters) <- c("altitude", "forestcover")
  coords <- terra::unwrap(melip.coords)
  surface <- conductance_surface(scale_covariates(rasters), coords,
                                  directions = 4, saveStack = TRUE)
  control <- NewtonRaphsonControl(maxit = 200)
  fit <- terradish(melip.Fst ~s(altitude, df = 3), surface,
                    smooth_loglinear_conductance, mlpe, control = control)
  coarse <- terradish(melip.Fst ~s(altitude, df = 3), surface,
                       smooth_loglinear_conductance, mlpe, control = control,
                       approximation = "coarse_raster",
                       approximation_control = list(factor = 2, exact_refine = TRUE))
  expect_equal(coef(fit), coef(coarse), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(coarse)), tolerance = 1e-8)
  clamped <- .clamp_graph_covariates(surface, support = "focal",
                                     support_probs = c(0, 1), clamp_covariates = NULL)
  unchanged <- surface$x$altitude == clamped$x$altitude
  rebuilt <- fit$submodels$f_factory(fit$formula, clamped$x)
  expect_equal(rebuilt(coef(fit))$conductance[unchanged],
                 fit$submodels$f(coef(fit))$conductance[unchanged], tolerance = 1e-10)
  none <- conductance(surface, fit, support = "none")
  focal <- conductance(surface, fit, support = "focal")
  values_none <- terra::values(none)[, 1]
  values_focal <- terra::values(focal)[, 1]
  expect_equal(values_focal[is.finite(values_focal)][unchanged],
                 values_none[is.finite(values_none)][unchanged], tolerance = 1e-10)
})
