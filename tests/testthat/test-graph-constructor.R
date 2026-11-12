test_that("the reduced graph Laplacian matches the sparse projection", {
  r <- terra::rast(nrows = 12, ncols = 12, nlyrs = 2,
                   xmin = 0, xmax = 12, ymin = 0, ymax = 12,
                   crs = "EPSG:3857")
  terra::values(r) <- cbind(seq_len(144), rep(1:12, each = 12))
  names(r) <- c("altitude", "forestcover")
  sites <- terra::xyFromCell(r, c(15, 34, 76, 118))
  g <- conductance_surface(r, sites, directions = 8)

  # Ground the last vertex in two algebraically equivalent ways. Keep the
  # former projection here as an independent reference for the constructor.
  Q <- g$laplacian
  N <- nrow(Q)
  full <- Q + Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In <- Matrix::Diagonal(N)[-N, , drop = FALSE]
  old <- Matrix::forceSymmetric(In %*% full %*% t(In))
  principal <- Matrix::forceSymmetric(full[-N, -N, drop = FALSE])
  expect_s4_class(principal, "symmetricMatrix")
  expect_s4_class(principal, "sparseMatrix")
  expect_equal(as.matrix(principal), as.matrix(old), tolerance = 1e-12)
  expect_identical(g$reduced_index, seq_len(N - 1L))

  # Recover the matrix represented by the constructor's permuted LDL factor.
  parts <- Matrix::expand(g$choleski)
  reconstructed <- crossprod(parts$P, tcrossprod(parts$L) %*% parts$P)
  expect_equal(as.matrix(reconstructed), as.matrix(old), tolerance = 1e-12)
  expect_equal(as.matrix(reconstructed), as.matrix(principal), tolerance = 1e-12)

  set.seed(111)
  rhs <- rnorm(N - 1L)
  reference <- Matrix::solve(Matrix::Cholesky(old, LDL = TRUE), rhs)
  expect_equal(as.matrix(Matrix::solve(g$choleski, rhs)),
               as.matrix(reference), tolerance = 1e-12)
  expect_equal(as.matrix(Matrix::solve(Matrix::Cholesky(principal, LDL = TRUE), rhs)),
               as.matrix(reference), tolerance = 1e-12)
})
