# vignettes/precompute.R
#
# Pre-computes slow simulation results for simulation-design.Rmd.
# Run this script manually whenever you want to refresh the artifacts:
#
#   source("vignettes/precompute.R")
#
# Output files (committed to source):
#   vignettes/vignette-nu-power.rds
#   vignettes/vignette-power.rds
#
# Increase nsim / n_designs for publication-quality estimates.

library(terradish)
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

surface <- conductance_surface(covariates, melip.coords,
                               directions = 8, saveStack = TRUE)

theta_true <- c(forestcover = 0.8, altitude = -0.5)
tau_true   <- 1.0
sigma_true <- 0.3

# ---- nu-power artifact --------------------------------------------------- #

detect_rate <- function(nu, nsim = 20, seed = 100) {
  set.seed(seed)
  hits <- logical(nsim)
  for (i in seq_len(nsim)) {
    s <- simulate_covariance_response(
      theta = theta_true, formula = ~ forestcover + altitude,
      data = surface, tau = tau_true, sigma = sigma_true, nu = nu
    )
    f <- terradish(
      s$covariance ~ forestcover + altitude, data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = wishart_covariance, nu = nu,
      control = NewtonRaphsonControl(maxit = 10, verbose = FALSE)
    )
    est <- coef(f)["altitude"]
    se  <- sqrt(diag(solve(f$fit$hessian)))["altitude"]
    hits[i] <- abs(est / se) > 1.96
  }
  mean(hits)
}

nu_grid     <- c(50, 100, 250, 500, 1000)
power_by_nu <- sapply(nu_grid, detect_rate)

saveRDS(
  list(nu_grid = nu_grid, power_by_nu = power_by_nu),
  "vignettes/vignette-nu-power.rds"
)
cat("Saved vignettes/vignette-nu-power.rds\n")

# ---- power-run artifact -------------------------------------------------- #

power <- covariance_response_power(
  theta             = theta_true,
  formula           = ~ forestcover + altitude,
  data              = surface,
  sample_sizes      = c(10, 20, 34),
  strategies        = c("spacefill", "random"),
  conductance_model = loglinear_conductance,
  fit_models = list(
    full          = list(formula = ~ forestcover + altitude),
    altitude_only = list(formula = ~ altitude)
  ),
  tau        = tau_true,
  sigma      = sigma_true,
  nu         = 500,
  nsim       = 10,
  n_designs  = 2,
  seed       = 1,
  control    = NewtonRaphsonControl(maxit = 10, verbose = FALSE)
)

saveRDS(power, "vignettes/vignette-power.rds")
cat("Saved vignettes/vignette-power.rds\n")
