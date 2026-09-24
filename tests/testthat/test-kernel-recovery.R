test_that("environmental kernel weights recover simulated covariance", {
  set.seed(251)
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 0, ymax = 8)
  terra::values(r) <- rnorm(terra::ncell(r))
  names(r) <- "x"
  coords <- terra::xyFromCell(r, c(1, 8, 13, 20, 27, 35, 42, 49, 57, 64))
  surface <- conductance_surface(r, coords, directions = 4)
  E <- terradish_distance(matrix(0.4, 1), ~x, surface, covariance = TRUE)$covariance[, , 1]
  g <- wishart_covariates(rnorm(nrow(E)), scale = TRUE)
  K <- attr(g, "kernel_covariates")[, , 1]
  Sigma <- 0.8 * E + 0.6 * K + 0.2 * diag(nrow(E))
  for (seed in 51:53) {
    set.seed(seed)
    S <- rWishart(1, 300, Sigma)[, , 1] / 300
    result <- radish_subproblem(g, E, S, nu = 300,
                                 control = NewtonRaphsonControl(maxit = 200))
    se <- sqrt(diag(solve(result$fit$hessian)))[2]
    expect_lt(abs(as.numeric(result$phi[2]) - 0.6), 3 * se)
    expect_lt(max(abs(result$fit$gradient)), 1e-4)
  }
})
