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

theta_true <- c(forestcover = 0.6, altitude = -0.4)
tau_true   <- 2
sigma_true <- 1

# ---- planning-screen artifact ------------------------------------------- #

altitude_surface <- conductance_surface(
  covariates[["altitude"]], melip.coords,
  directions = 8, saveStack = TRUE
)

screen_theta   <- c(altitude = -0.25)
screen_sigma   <- 1
signal_ratios  <- c(0.5, 1, 2, 4)
screen_nsim    <- 10

screen_power <- function(sample_size, nu, signal_ratio, seed) {
  assessment <- covariance_response_power(
    theta             = screen_theta,
    formula           = ~ altitude,
    data              = altitude_surface,
    sample_sizes      = sample_size,
    strategies        = "spacefill",
    tau               = signal_ratio * screen_sigma,
    sigma             = screen_sigma,
    nu                = nu,
    nsim              = screen_nsim,
    seed              = seed,
    control           = NewtonRaphsonControl(maxit = 10, verbose = FALSE)
  )
  out <- assessment$parameter_summary
  stopifnot(nrow(out) == 1L, out$parameter == "altitude")
  data.frame(
    sample_size = sample_size,
    nu = nu,
    signal_ratio = signal_ratio,
    fit_rate = out$fit_rate,
    power = out$power
  )
}

screen_grid <- function(sample_sizes, nu_values, seed) {
  grid <- expand.grid(
    sample_size = sample_sizes,
    nu = nu_values,
    signal_ratio = signal_ratios
  )
  do.call(
    rbind,
    Map(
      screen_power,
      sample_size = grid$sample_size,
      nu = grid$nu,
      signal_ratio = grid$signal_ratio,
      seed = seed + seq_len(nrow(grid))
    )
  )
}

marker_power <- screen_grid(sample_sizes = 20, nu_values = c(20, 50, 100, 200),
                            seed = 100)
site_power <- screen_grid(sample_sizes = c(8, 12, 20, 34), nu_values = 100,
                          seed = 200)
saveRDS(
  list(
    marker_power = marker_power,
    site_power = site_power,
    signal_ratios = signal_ratios,
    screen_nsim = screen_nsim
  ),
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
  nu         = 100,
  nsim       = 20,
  n_designs  = 2,
  seed       = 1,
  control    = NewtonRaphsonControl(maxit = 10, verbose = FALSE)
)

saveRDS(power, "vignettes/vignette-power.rds")
cat("Saved vignettes/vignette-power.rds\n")
