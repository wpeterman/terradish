# Generate the precomputed results object for the hierarchical-conductance vignette.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))
set.seed(20)

data(melip)
covs <- c(terra::scale(terra::unwrap(melip.altitude)),
          terra::scale(terra::unwrap(melip.forestcover)))
names(covs) <- c("altitude", "forestcover")
surface <- conductance_surface(covs, terra::unwrap(melip.coords), directions = 8)
n  <- length(surface$demes)
Vc <- surface$vertex_coordinates
X  <- as.matrix(surface$x[, c("altitude", "forestcover")])

xr <- range(Vc[, 1]); yr <- range(Vc[, 2])
blob <- 1.6 * exp(-(((Vc[, 1] - (xr[1] + 0.35 * diff(xr)))^2 +
                     (Vc[, 2] - (yr[1] + 0.6 * diff(yr)))^2) /
                    (2 * (0.18 * sqrt(diff(xr)^2 + diff(yr)^2))^2)))
blob <- blob - mean(blob)
E <- as.matrix(terradish_algorithm(
  terradish:::.design_loglinear_model(cbind(X, blob)), leastsquares, surface,
  S = diag(n), theta = c(0.5, -0.6, 1.0),
  objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
nu <- 1500
S <- rWishart(1, df = nu, Sigma = (E + 0.25 * diag(n)) / nu)[, , 1]

fit_h <- terradish_hierarchical(S ~ altitude + forestcover, data = surface,
                                measurement_model = wishart_covariance, nu = nu,
                                field_resolution = 6, tau2 = "reml", verbose = FALSE)
fit_plain <- terradish(S ~ altitude + forestcover, data = surface,
                       conductance_model = loglinear_conductance,
                       measurement_model = wishart_covariance, nu = nu, leverage = FALSE,
                       control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))

coef_table <- rbind(
  `hierarchical (with field)` = c(coef(fit_h)),
  `plain terradish (no field)` = c(coef(fit_plain)))
colnames(coef_table) <- names(coef(fit_h))

u_cells <- as.vector(fit_h$field$Z %*% fit_h$u)
field_raster <- terra::wrap(conductance_field(fit_h, surface, type = "field"))

res <- list(
  coef_table = coef_table,
  cor_field_blob = cor(u_cells, blob),
  tau2 = fit_h$tau2,
  tau2_selection = fit_h$tau2_selection,
  field_raster = field_raster)

saveRDS(res, "vignettes/vignette-hierarchical.rds")
cat("saved vignettes/vignette-hierarchical.rds\n")
cat(sprintf("cor(field, blob) = %.3f; tau2 = %.4g\n", res$cor_field_blob, res$tau2))
print(round(coef_table, 3))
