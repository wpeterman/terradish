brute_window_gaussian <- function(mat, sigma) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  out <- mat
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    if (!is.finite(mat[i, j])) next
    numerator <- denominator <- 0
    for (a in seq_len(nr)) for (b in seq_len(nc)) {
      kr <- i - a + ceiling(nr / 2)
      kc <- j - b + ceiling(nc / 2)
      if (kr < 1 || kr > nr || kc < 1 || kc > nc || !is.finite(mat[a, b])) next
      weight <- exp(-((i - a)^2 + (j - b)^2) / (2 * sigma^2))
      numerator <- numerator + weight * mat[a, b]
      denominator <- denominator + weight
    }
    out[i, j] <- numerator / denominator
  }
  out
}

test_that("Gaussian convolution is aligned on odd and even rasters", {
  set.seed(401)
  for (dims in list(c(7, 9), c(8, 10))) {
    r <- terra::rast(nrows = dims[1], ncols = dims[2], xmin = 0,
                     xmax = dims[2], ymin = 0, ymax = dims[1])
    terra::values(r) <- rnorm(terra::ncell(r))
    r[c(1, 14)] <- NA
    rc <- terra::rowColFromCell(r, seq_len(terra::ncell(r)))
    prep <- .gaussian_scale_prepare_layer(r, rc)
    mat <- as.matrix(r, wide = TRUE)
    for (sigma in c(0.7, 1.5, 3)) {
      result <- .gaussian_scale_layer_values(prep, sigma, standardize = FALSE)
      truth <- brute_window_gaussian(mat, sigma)
      expect_equal(result$value, truth[rc], tolerance = 1e-10)
      keep <- is.finite(result$value)
      numerical <- numDeriv::jacobian(function(s)
        .gaussian_scale_layer_values(prep, s, standardize = FALSE)$value[keep], sigma)
      expect_equal(result$deriv[keep], c(numerical), tolerance = 1e-7)
    }
  }
  terra::values(r) <- 0
  r[terra::cellFromRowCol(r, 4, 5)] <- 1
  spike <- .gaussian_scale_layer_values(.gaussian_scale_prepare_layer(r, rc),
                                        0.7, standardize = FALSE)$value
  expect_equal(which.max(spike), terra::cellFromRowCol(r, 4, 5))
})

test_that("Gaussian upper bounds respect the retained kernel window", {
  r <- terra::rast(nrows = 12, ncols = 18, xmin = 0, xmax = 1800,
                   ymin = 0, ymax = 1200)
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "x"
  bounds <- .gaussian_scale_sigma_bounds(list(stack = r), "x")
  expect_equal(unname(bounds$upper), 200)
  expect_warning(.gaussian_scale_sigma_bounds(list(stack = r), "x", sigma_upper = 201),
                  "three-sigma support")
})
