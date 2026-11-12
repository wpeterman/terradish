# The two melip rasters used to ship with one shared layer name. They are now
# named for what they hold, and nothing else about the data changed.

test_that("melip rasters carry distinct layer names", {
  data(melip, package = "terradish", envir = environment())
  altitude <- terra::unwrap(melip.altitude)
  forestcover <- terra::unwrap(melip.forestcover)

  expect_identical(names(altitude), "altitude")
  expect_identical(names(forestcover), "forestcover")
  expect_identical(names(c(altitude, forestcover)),
                   c("altitude", "forestcover"))
})

test_that("melip values, grid, sites, and FST are as documented", {
  data(melip, package = "terradish", envir = environment())
  altitude <- terra::unwrap(melip.altitude)
  forestcover <- terra::unwrap(melip.forestcover)
  coords <- terra::unwrap(melip.coords)

  expect_identical(dim(altitude), c(156, 180, 1))
  expect_identical(dim(forestcover), c(156, 180, 1))
  expect_true(terra::compareGeom(altitude, forestcover))
  expect_equal(as.vector(terra::ext(altitude)),
               c(-42.9125, -34.8125, -9.4975, -2.4775),
               ignore_attr = TRUE)
  expect_true(terra::is.lonlat(altitude))
  expect_identical(sum(!is.na(terra::values(altitude))), 22443L)
  expect_equal(range(terra::values(forestcover), na.rm = TRUE),
               c(0.004, 0.948), tolerance = 0.01)

  expect_identical(nrow(coords), 37)
  expect_identical(dim(melip.Fst), c(37L, 37L))
  expect_true(isSymmetric(melip.Fst))
  expect_equal(range(melip.Fst), c(-0.062, 0.280), tolerance = 0.01)
})

test_that("renaming the layers, as older scripts do, still works", {
  data(melip, package = "terradish", envir = environment())
  covariates <- c(terra::unwrap(melip.altitude),
                  terra::unwrap(melip.forestcover))
  names(covariates) <- c("altitude", "forestcover")
  scaled <- scale_covariates(covariates)
  expect_named(attr(scaled, "terradish_scale"), c("altitude", "forestcover"))

  # Without renaming, the stored scaling is keyed by the shipped names.
  unrenamed <- scale_covariates(c(terra::unwrap(melip.altitude),
                                  terra::unwrap(melip.forestcover)))
  expect_named(attr(unrenamed, "terradish_scale"),
               c("altitude", "forestcover"))
})
