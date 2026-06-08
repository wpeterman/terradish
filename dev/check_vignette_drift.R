# Smoke-test the wishart-covariance vignette drift-surface chunks end to end.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)
covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)
surface <- conductance_surface(covariates, melip.coords, directions = 8)

# --- drift-simulate chunk ---
n_site <- length(surface$demes)
site_drift <- scale(seq_len(n_site))[, 1]
E_true <- as.matrix(terradish_algorithm(
  loglinear_conductance(~ forestcover + altitude, surface$x),
  leastsquares, surface, S = diag(n_site),
  theta = c(forestcover = 1.0, altitude = -0.5),
  objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE
)$covariance)
Z <- cbind(1, scale(site_drift, scale = FALSE))
gamma_true <- c(log(0.2), -0.6)
nugget     <- as.vector(exp(Z %*% gamma_true))
Sigma_true <- 1.0 * E_true + diag(nugget)
set.seed(123)
nu_demo <- 500
S_drift <- rWishart(1, df = nu_demo, Sigma = Sigma_true / nu_demo)[, , 1]

# --- drift-fit chunk ---
g_drift <- wishart_drift_covariates(site_drift, model = "wishart_covariance")
fit_drift <- suppressWarnings(terradish(S_drift ~ forestcover + altitude, data = surface,
  conductance_model = loglinear_conductance, measurement_model = g_drift, nu = nu_demo))
fit_scalar <- suppressWarnings(terradish(S_drift ~ forestcover + altitude, data = surface,
  conductance_model = loglinear_conductance, measurement_model = wishart_covariance, nu = nu_demo))

cat("phi rownames:", rownames(fit_drift$fit$phi), "\n")
cat("gamma_var1 (true -0.6):", round(fit_drift$fit$phi["gamma_var1", 1], 3), "\n")
aic <- function(f) -2 * f$loglik + 2 * (length(coef(f)) + nrow(f$fit$phi))
cat("AIC drift:", round(aic(fit_drift), 1), " scalar:", round(aic(fit_scalar), 1),
    " drift better:", aic(fit_drift) < aic(fit_scalar), "\n")
cat("summary() works:\n")
print(summary(fit_drift))
cat("\nVIGNETTE DRIFT CHUNKS OK\n")
