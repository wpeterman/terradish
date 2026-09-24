.libPaths(c("dev/check/phase2-library", .libPaths()))
pkgload::load_all("dev/check/phase5-package", quiet = TRUE)
source("dev/check/phase5-package/tests/testthat/helper-robustness.R")
g <- robustness_surface()
g <- conductance_surface(g$stack, terra::xyFromCell(g$stack,
  c(1, 5, 8, 12, 20, 27, 33, 37, 45, 53, 57, 64)))
coords <- g$vertex_coordinates[g$demes, ]
z <- pairwise_covariates(pairwise_endpoint_covariates(data.frame(climate = sin(1:12))),
                        geographic = dist(coords))
S <- simulate_covariance_response(c(x = 0.8), ~x, g, nu = 1000, seed = 301)$S
kernels <- attr(wishart_covariates(z, model = "wishart_covariance"), "kernel_covariates")
S <- S + 0.1 * kernels[, , 1] + 0.05 * kernels[, , 2]
D <- dist_from_cov(S)
fits <- list()
for (optimizer in c("newton", "bfgs")) {
  fits[[optimizer]] <- terradish(D ~ x, g, measurement_model = mlpe_covariates(z),
    optimizer = optimizer, control = NewtonRaphsonControl(maxit = 100))
  print(fits[[optimizer]]$convergence)
  print(fits[[optimizer]]$mle$theta)
  print(fits[[optimizer]]$fit$phi)
}
saveRDS(fits, "dev/release_0.1.0/phase5_pairwise_diagnosis.rds")
