test_that("conductance_surface and conductance work with terra inputs", {
  fx <- fit_fixture(control = NewtonRaphsonControl(maxit = 2, verbose = FALSE))

  expect_s3_class(fx$surface, "terradish_graph")
  expect_s3_class(fx$surface, "radish_graph")
  expect_true(inherits(fx$surface$stack, "SpatRaster"))

  cond <- conductance(fx$surface, fx$fit)
  expect_true(inherits(cond, "SpatRaster"))
  expect_equal(terra::nlyr(cond), 3)
})

test_that("utility helpers behave as expected", {
  expect_equal(scale_to_0_1(c(2, 4, 6)), c(0, 0.5, 1))
  expect_equal(lower(matrix(1:9, nrow = 3)), c(2L, 3L, 6L))

  r <- terra::rast(nrows = 2, ncols = 2, vals = c(1, 3, 5, NA))
  out <- scale_to_0_1(r)
  expect_true(inherits(out, "SpatRaster"))
  expect_equal(as.vector(terra::values(out)), c(0, 0.5, 1, NA))

  covs <- c(r, 2 * r + 1)
  names(covs) <- c("a", "b")
  scaled_covs <- scale_covariates(covs)
  expect_true(inherits(scaled_covs, "SpatRaster"))
  expect_equal(names(attr(scaled_covs, "terradish_scale")), c("a", "b"))
  expect_equal(attr(scaled_covs, "terradish_scale")$a[["center"]], 3)
  expect_equal(attr(scaled_covs, "terradish_scale")$a[["scale"]], 2)

  ranged_covs <- scale_covariates(covs, method = "minmax")
  expect_equal(as.vector(terra::values(ranged_covs[[1]])), c(0, 0.5, 1, NA))
  expect_equal(attr(ranged_covs, "terradish_scale")$a[["center"]], 1)
  expect_equal(attr(ranged_covs, "terradish_scale")$a[["scale"]], 4)
})

test_that("crop_to_focal_buffer reduces raster extent before graph construction", {
  r <- terra::rast(nrows = 20, ncols = 20, vals = seq_len(400),
                   xmin = 0, xmax = 20, ymin = 0, ymax = 20)
  covs <- c(r, r * 2)
  names(covs) <- c("a", "b")
  covs <- scale_covariates(covs)
  pts <- matrix(c(4.5, 4.5,
                  6.5, 6.5,
                  7.5, 5.5),
                ncol = 2, byrow = TRUE)

  cropped <- crop_to_focal_buffer(covs, pts, buffer = 2)
  expect_s4_class(cropped, "SpatRaster")
  expect_lt(terra::ncell(cropped), terra::ncell(covs))
  expect_equal(names(attr(cropped, "terradish_scale")), c("a", "b"))
  expect_false(anyNA(terra::cellFromXY(cropped[[1]], pts)))

  full_surface <- conductance_surface(covs, pts, directions = 4)
  cropped_surface <- conductance_surface(covs, pts, directions = 4,
                                         crop_buffer = 2)
  expect_lt(nrow(cropped_surface$x), nrow(full_surface$x))
  expect_equal(length(cropped_surface$demes), nrow(pts))
})

test_that("pca_dist works when adegenet is available", {
  skip_if_not_installed("adegenet")
  data(nancycats, package = "adegenet")
  out <- pca_dist(nancycats, n_axes = 4)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(adegenet::nInd(nancycats), adegenet::nInd(nancycats)))
})
