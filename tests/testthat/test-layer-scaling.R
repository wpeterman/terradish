test_that("unit scaling preserves each raster layer's own range", {
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = 3)
  terra::values(r) <- cbind(c(1, 3, 5, NA), c(100, 150, 200, NA), c(7, 7, 7, NA))
  out <- scale_to_0_1(r)
  expect_equal(terra::values(out),
    cbind(c(0, 0.5, 1, NA), c(0, 0.5, 1, NA), c(0, 0, 0, NA)),
    ignore_attr = TRUE)
})
