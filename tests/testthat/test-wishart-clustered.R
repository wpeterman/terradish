test_that("contrast Wishart recovers clustered-site simulations", {
  set.seed(2)
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30)
  f1 <- f2 <- r
  terra::values(f1) <- rnorm(terra::ncell(r))
  f1 <- terra::focal(f1, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  terra::values(f2) <- rnorm(terra::ncell(r))
  f2 <- terra::focal(f2, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  covariates <- c(f1, f2)
  names(covariates) <- c("x1", "x2")
  covariates <- scale_covariates(covariates)
  xy <- cbind(runif(20, 1, 12), runif(20, 1, 12))
  surface <- suppressWarnings(conductance_surface(covariates, xy, directions = 8))
  truth <- c(x1 = 0.6, x2 = -0.4)
  n <- length(surface$demes)
  H <- diag(n) - matrix(1 / n, n, n)
  L <- qr.Q(qr(stats::contr.helmert(n)))
  for (seed in 101:103) {
    sim <- simulate_covariance_response(truth, ~x1 + x2, surface,
                                         tau = 2, sigma = 0.5, nu = 300, seed = seed)
    S <- H %*% sim$covariance %*% H
    D <- outer(diag(S), diag(S), "+") - 2 * S
    control <- NewtonRaphsonControl(maxit = 200)
    a <- terradish(S ~x1 + x2, surface, loglinear_conductance,
                   wishart_covariance, nu = 300, control = control)
    b <- terradish(D ~x1 + x2, surface, loglinear_conductance,
                   generalized_wishart, nu = 300, control = control)
    expect_gt(as.numeric(a$fit$phi[1]), 0)
    se <- sqrt(diag(-solve(a$mle$hessian)))
    expect_true(all(abs(coef(a) - truth) < 3 * se))
    expect_equal(coef(a), coef(b), tolerance = 1e-6)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-6)
    expect_equal(crossprod(L, a$fit$fitted %*% L),
                  -0.5 * crossprod(L, b$fit$fitted %*% L), tolerance = 1e-6)
  }
})
