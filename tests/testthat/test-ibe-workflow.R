test_that("pairwise endpoint covariates can be built from raster values", {
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.coords <- terra::unwrap(melip.coords)

  names(melip.altitude) <- "altitude"
  keep <- 1:10
  z <- expect_no_warning(
    pairwise_endpoint_covariates(melip.altitude,
                                 melip.coords[keep],
                                 transform = "absdiff",
                                 scale = TRUE)
  )

  expect_s3_class(z, "radish_pairwise_covariates")
  expect_equal(dim(z), c(choose(length(keep), 2), 1))
  expect_equal(colnames(z), "absdiff_altitude")
})

test_that("joint IBE + IBR MLPE models fit through radish", {
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  keep <- 1:10
  melip.Fst_small <- melip.Fst[keep, keep, drop = FALSE]
  coords <- melip.coords[keep]

  names(melip.altitude) <- "altitude_site"
  covariates <- c(terra::scale(melip.altitude),
                  terra::scale(melip.forestcover))
  names(covariates) <- c("altitude", "forestcover")
  surface <- conductance_surface(covariates, coords, directions = 8)

  g_joint <- mlpe_covariates(melip.altitude, coords,
                             transform = "absdiff",
                             scale = TRUE)

  fit <- expect_no_warning(
    radish(melip.Fst_small ~ altitude + forestcover,
           data = surface,
           conductance_model = loglinear_conductance,
           measurement_model = g_joint,
           leverage = FALSE,
           control = NewtonRaphsonControl(maxit = 8, verbose = FALSE))
  )

  expect_s3_class(fit, "radish")
  expect_true(is.finite(fit$loglik))

  sm <- summary(fit)
  expect_s3_class(sm, "summary.radish")
  expect_true(any(grepl("absdiff_altitude_site", names(sm$phi))))
})

test_that("radish_cv rebuilds MLPE endpoint-covariate models on splits", {
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  keep <- 1:12
  melip.Fst_small <- melip.Fst[keep, keep, drop = FALSE]
  coords <- melip.coords[keep]

  names(melip.altitude) <- "altitude_site"
  covariates <- c(terra::scale(melip.altitude),
                  terra::scale(melip.forestcover))
  names(covariates) <- c("altitude", "forestcover")

  g_joint <- mlpe_covariates(melip.altitude, coords,
                             transform = "absdiff",
                             scale = TRUE)

  out <- expect_no_warning(
    radish_cv(coords,
              covariates,
              melip.Fst_small ~ altitude + forestcover,
              model = g_joint,
              prop_train = 2 / 3,
              seed = 1,
              fit_full = FALSE,
              control = NewtonRaphsonControl(maxit = 8, verbose = FALSE))
  )

  expect_true(is.list(out))
  expect_s3_class(out$train_mod, "radish")
  expect_length(out$cv_loglik, 1L)
  expect_true(is.finite(out$cv_loglik))
})
