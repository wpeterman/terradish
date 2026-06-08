# Exercise the hardened terradish_hierarchical() API end-to-end.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))
set.seed(20)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)
covs <- c(terra::scale(melip.altitude), terra::scale(melip.forestcover))
names(covs) <- c("altitude", "forestcover")
surface <- conductance_surface(covs, melip.coords, directions = 8)
n  <- length(surface$demes)
Vc <- surface$vertex_coordinates
X  <- as.matrix(surface$x[, c("altitude", "forestcover")])

# TRUE surface: mapped covariates + UNMAPPED smooth blob
loglinear_from_matrix <- terradish:::.design_loglinear_model
xr <- range(Vc[, 1]); yr <- range(Vc[, 2])
bx <- xr[1] + 0.35 * diff(xr); by <- yr[1] + 0.6 * diff(yr)
rad <- 0.18 * sqrt(diff(xr)^2 + diff(yr)^2)
blob <- 1.6 * exp(-(((Vc[, 1] - bx)^2 + (Vc[, 2] - by)^2) / (2 * rad^2)))
blob <- blob - mean(blob)
theta_true <- c(altitude = 0.5, forestcover = -0.6)
E_true <- as.matrix(terradish_algorithm(
  loglinear_from_matrix(cbind(X, blob)), leastsquares, surface, S = diag(n),
  theta = c(theta_true, 1.0), objective = FALSE,
  gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
nu <- 1500
S <- rWishart(1, df = nu, Sigma = (1.0 * E_true + 0.25 * diag(n)) / nu)[, , 1]

cat("=== (A) terradish_hierarchical (tau2 fixed at REML optimum 0.3162) ===\n")
fit <- terradish_hierarchical(S ~ altitude + forestcover, data = surface,
                              measurement_model = wishart_covariance, nu = nu,
                              field_resolution = 6, tau2 = 0.3162, verbose = FALSE)
print(fit)
u_cells <- as.vector(fit$field$Z %*% fit$u)
cat(sprintf("\nselected tau2 = %.4g\n", fit$tau2))
cat(sprintf("cor(field, true blob) = %.3f\n", cor(u_cells, blob)))
cat(sprintf("theta: alt=%.3f (0.50), fc=%.3f (-0.60)\n", fit$theta[1], fit$theta[2]))

cat("\n=== (B) no-field terradish baseline (omitted-variable bias) ===\n")
fit_nf <- terradish(S ~ altitude + forestcover, data = surface,
                    conductance_model = loglinear_conductance,
                    measurement_model = wishart_covariance, nu = nu, leverage = FALSE,
                    control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))
cat(sprintf("no-field theta: alt=%.3f, fc=%.3f\n", coef(fit_nf)[1], coef(fit_nf)[2]))

cat("\n=== (C) conductance_field() returns a raster ===\n")
fld_r <- conductance_field(fit, surface, type = "field")
cat("class:", class(fld_r)[1], " nlyr:", terra::nlyr(fld_r), " name:", names(fld_r), "\n")
logc_r <- conductance_field(fit, surface, type = "logconductance")
cat("logconductance raster OK:", inherits(logc_r, "SpatRaster"), "\n")

cat("\n=== (D) reduction: tiny tau2 shrinks the field ~ terradish ===\n")
fit_small <- terradish_hierarchical(S ~ altitude + forestcover, data = surface,
                                    measurement_model = wishart_covariance, nu = nu,
                                    field_resolution = 6, tau2 = 1e-4, verbose = FALSE)
cat(sprintf("tau2=1e-4: max|u| = %.4f; theta alt=%.3f fc=%.3f (vs no-field %.3f %.3f)\n",
            max(abs(fit_small$u)), fit_small$theta[1], fit_small$theta[2],
            coef(fit_nf)[1], coef(fit_nf)[2]))

cat("\n=== (E) summary() works ===\n")
invisible(summary(fit))

cat("\nDONE\n")
