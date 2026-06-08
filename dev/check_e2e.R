suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))
set.seed(7)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(terra::scale(melip.altitude), terra::scale(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")

## ---- (B) gradient validation on a small graph (fast finite differences) ----
keepB <- 1:10
surfB <- conductance_surface(covariates, melip.coords[keepB, ], directions = 8)
nB <- length(surfB$demes)
cm_B <- loglinear_conductance(~ altitude + forestcover, surfB$x)
ZsiteB <- scale(seq_len(nB))[, 1]
gB <- wishart_drift_covariates(ZsiteB, model = "wishart_covariance")

theta_true <- c(altitude = 0.2, forestcover = -0.25)
alg0 <- terradish_algorithm(cm_B, leastsquares, surfB, S = diag(nB),
                            theta = theta_true, objective = FALSE,
                            gradient = FALSE, hessian = FALSE, partial = FALSE)
E_true <- as.matrix(alg0$covariance)
Zd <- cbind(1, scale(ZsiteB, center = TRUE, scale = FALSE))
nug <- as.vector(exp(Zd %*% c(log(0.2), 0.6)))
nu <- 200
SB <- rWishart(1, df = nu, Sigma = (0.8 * E_true + diag(nug)) / nu)[, , 1]

theta0 <- c(altitude = 0.1, forestcover = -0.15)
prof_obj <- function(th)
  terradish_algorithm(cm_B, gB, surfB, S = SB, theta = th, nu = nu,
                      gradient = FALSE, hessian = FALSE, partial = FALSE)$objective
analytic <- terradish_algorithm(cm_B, gB, surfB, S = SB, theta = theta0, nu = nu,
                                hessian = FALSE, partial = FALSE)$gradient
h <- 1e-4
num <- sapply(seq_along(theta0), function(k) {
  tp <- theta0; tm <- theta0; tp[k] <- tp[k] + h; tm[k] <- tm[k] - h
  (prof_obj(tp) - prof_obj(tm)) / (2 * h)
})
cat("=== (B) profiled theta-gradient: analytic vs finite-difference ===\n")
cat("analytic:", round(analytic, 5), "\n")
cat("numeric :", round(num, 5), "\n")
cat("max abs diff:", max(abs(analytic - num)),
    ifelse(max(abs(analytic - num)) < 1e-4, " OK", " **FAIL**"), "\n")

## ---- (C) recovery + model selection on a larger graph ----
keepC <- 1:28
surfC <- conductance_surface(covariates, melip.coords[keepC, ], directions = 8)
nC <- length(surfC$demes)
cm_C <- loglinear_conductance(~ altitude + forestcover, surfC$x)
## drift covariate = scaled altitude AT the focal sites (a real site covariate)
site_alt <- terra::extract(covariates[["altitude"]], melip.coords[keepC, ])[, 2]
Zc <- scale(site_alt)[, 1]

alg0C <- terradish_algorithm(cm_C, leastsquares, surfC, S = diag(nC),
                             theta = theta_true, objective = FALSE,
                             gradient = FALSE, hessian = FALSE, partial = FALSE)
E_C <- as.matrix(alg0C$covariance)
ZdC <- cbind(1, scale(Zc, center = TRUE, scale = FALSE))
gamma_true <- c(log(0.25), 0.7)
nugC <- as.vector(exp(ZdC %*% gamma_true))
nuC <- 500
Sigma_C <- 0.8 * E_C + diag(nugC)

cat("\n=== (C) recovery of drift slope + AIC vs scalar nugget (5 draws) ===\n")
gam_hats <- numeric(0); win <- 0
for (rep in 1:5) {
  S_C <- rWishart(1, df = nuC, Sigma = Sigma_C / nuC)[, , 1]
  gC <- wishart_drift_covariates(Zc, model = "wishart_covariance")
  fit_d <- suppressWarnings(terradish(S_C ~ altitude + forestcover, data = surfC,
                     conductance_model = loglinear_conductance,
                     measurement_model = gC, nu = nuC, leverage = FALSE,
                     control = NewtonRaphsonControl(maxit = 100, verbose = FALSE)))
  fit_s <- suppressWarnings(terradish(S_C ~ altitude + forestcover, data = surfC,
                     conductance_model = loglinear_conductance,
                     measurement_model = wishart_covariance, nu = nuC, leverage = FALSE,
                     control = NewtonRaphsonControl(maxit = 100, verbose = FALSE)))
  kd <- length(fit_d$mle$theta) + nrow(fit_d$fit$phi)
  ks <- length(fit_s$mle$theta) + nrow(fit_s$fit$phi)
  aic_d <- -2 * fit_d$loglik + 2 * kd
  aic_s <- -2 * fit_s$loglik + 2 * ks
  gh <- fit_d$fit$phi["gamma_var1", 1]
  gam_hats <- c(gam_hats, gh)
  if (aic_d < aic_s) win <- win + 1
  cat(sprintf("  rep %d: gamma_hat=%.3f (true %.2f)  AIC drift=%.1f scalar=%.1f  drift_wins=%s\n",
              rep, gh, gamma_true[2], aic_d, aic_s, aic_d < aic_s))
}
cat(sprintf("mean gamma_hat=%.3f (true %.2f); drift preferred %d/5\n",
            mean(gam_hats), gamma_true[2], win))
cat("\nDONE\n")
