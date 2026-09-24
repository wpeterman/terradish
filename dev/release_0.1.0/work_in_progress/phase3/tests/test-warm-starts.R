test_that("nuisance vectors are used and named, with an audible fallback", {
  set.seed(17)
  X <- matrix(rnorm(25), 5)
  E <- crossprod(X)
  S <- .8 * E + diag(5)
  seen <- list()
  model <- function(E, S, phi, ...) {
    if (missing(phi)) return(wishart_covariance(E, S))
    seen[[length(seen) + 1L]] <<- phi
    if (phi[1] == 9999) stop("deliberately invalid warm start")
    wishart_covariance(E, S, phi = phi, ...)
  }
  warm <- c(.7, log(.9))
  result <- radish_subproblem(model, E, S, nu = 30, phi = warm,
    control = NewtonRaphsonControl(verbose = FALSE))
  expect_equal(unname(seen[[1]]), warm)
  expect_equal(rownames(result$phi), c("tau", "sigma"))
  expect_message(radish_subproblem(model, E, S, nu = 30, phi = c(9999, 0),
    control = NewtonRaphsonControl(verbose = TRUE)), "Nuisance warm start failed")
})

test_that("a poor nondefault start on the no-structure boundary is retried once", {
  set.seed(4)
  r <- terra::rast(nrows = 25, ncols = 25, xmin = 0, xmax = 25, ymin = 0, ymax = 25)
  terra::values(r) <- rnorm(625)
  x1 <- terra::focal(r, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  terra::values(r) <- rnorm(625)
  x2 <- terra::focal(r, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  rasters <- c(x1, x2); names(rasters) <- c("x1", "x2")
  rasters <- scale_covariates(rasters)
  coords <- cbind(runif(15, 1, 24), runif(15, 1, 24))
  graph <- conductance_surface(rasters, coords, directions = 8)
  sim <- simulate_covariance_response(c(x1 = .7, x2 = -.3), ~x1 + x2,
    graph, tau = 1, sigma = .2, nu = 200, seed = 1)
  D <- dist_from_cov(sim$covariance)
  fit <- terradish(D ~ x1 + x2, graph, measurement_model = mlpe, theta = c(2, 2))
  expect_true(fit$convergence$restarted)
  expect_false(.no_structure_boundary(fit$fit))
  expect_gt(fit$loglik, 260)
})
