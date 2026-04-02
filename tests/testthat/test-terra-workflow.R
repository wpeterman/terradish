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
})

test_that("pca_dist works when adegenet is available", {
  skip_if_not_installed("adegenet")
  data(nancycats, package = "adegenet")
  out <- pca_dist(nancycats, n_axes = 4)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(adegenet::nInd(nancycats), adegenet::nInd(nancycats)))
})
