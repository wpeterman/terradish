test_that("Gaussian post-smoothing standardization cancels affine pre-scaling", {
  layer <- terra::rast(nrows = 5, ncols = 5, vals = seq(2, 50, by = 2))
  layer[7] <- NA_real_
  names(layer) <- "covariate"

  cells <- which(is.finite(terra::values(layer)[, 1]))
  rowcol <- terra::rowColFromCell(layer, cells)
  original <- .gaussian_scale_prepare_layer(layer, rowcol)
  reference <- .gaussian_scale_layer_values(original, sigma = 1.5)

  for (method in c("zscore", "minmax")) {
    prescaled <- scale_covariates(layer, method = method)
    prepared <- .gaussian_scale_prepare_layer(prescaled, rowcol)
    result <- .gaussian_scale_layer_values(prepared, sigma = 1.5)

    expect_equal(result$value, reference$value, tolerance = 1e-12)
    expect_equal(result$deriv, reference$deriv, tolerance = 1e-12)
    expect_equal(result$second, reference$second, tolerance = 1e-12)
  }
})
