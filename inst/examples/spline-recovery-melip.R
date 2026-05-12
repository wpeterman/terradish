# Recoverability check for spline conductance models
#
# This script is intentionally small enough to run from RStudio. It has two
# parts:
# 1. Simulate a covariance response from a known non-linear conductance surface
#    and compare linear vs. spline recovery.
# 2. Fit the same first-pass spline model to the built-in melip data.
#
# After installation, find this file with:
# system.file("examples", "spline-recovery-melip.R", package = "terradish")

library(terradish)
library(terra)

run_spline_recovery_example <- function(keep = 1:10,
                                        nu = 150,
                                        maxit = 12,
                                        seed = 11)
{
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- scale_covariates(covariates)

  surface <- conductance_surface(covariates, melip.coords[keep], directions = 8)
  true_formula <- ~ forestcover + s(altitude, df = 4)
  true_model <- smooth_loglinear_conductance(true_formula, surface$x)
  theta_true <- attr(true_model, "default")
  theta_true["forestcover"] <- -0.15
  theta_true[grepl("^s\\(altitude\\)", names(theta_true))] <-
    c(-0.65, 1.1, -0.9, 0.45)

  sim <- simulate_covariance_response(
    theta = theta_true,
    formula = true_formula,
    data = surface,
    conductance_model = smooth_loglinear_conductance,
    tau = 0.8,
    sigma = 0.12,
    nu = nu,
    seed = seed
  )
  genetic_cov <- sim$covariance

  control <- NewtonRaphsonControl(maxit = maxit, verbose = FALSE)
  fit_linear <- suppressWarnings(
    terradish(
      genetic_cov ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = wishart_covariance,
      nu = nu,
      leverage = FALSE,
      control = control
    )
  )
  fit_spline <- suppressWarnings(
    terradish(
      genetic_cov ~ forestcover + s(altitude, df = 4),
      data = surface,
      conductance_model = smooth_loglinear_conductance,
      measurement_model = wishart_covariance,
      nu = nu,
      leverage = FALSE,
      control = control
    )
  )

  true_cond <- true_model(theta_true)$conductance
  linear_cond <- fit_linear$submodels$f(fit_linear$mle$theta)$conductance
  spline_cond <- fit_spline$submodels$f(fit_spline$mle$theta)$conductance

  model_comparison <- data.frame(
    model = c("linear", "spline"),
    logLik = c(as.numeric(logLik(fit_linear)), as.numeric(logLik(fit_spline))),
    AIC = c(AIC(fit_linear), AIC(fit_spline)),
    conductance_cor = c(cor(true_cond, linear_cond),
                        cor(true_cond, spline_cond)),
    conductance_spearman = c(cor(true_cond, linear_cond, method = "spearman"),
                             cor(true_cond, spline_cond, method = "spearman"))
  )

  list(
    truth = list(formula = true_formula, theta = theta_true,
                 conductance = true_cond),
    response = genetic_cov,
    surface = surface,
    fits = list(linear = fit_linear, spline = fit_spline),
    model_comparison = model_comparison,
    spline_conductance_raster = conductance(surface, fit_spline),
    spline_marginal_plot = plot(fit_spline, type = "marginal",
                                data = surface, n = 40)
  )
}

fit_melip_spline_example <- function(keep = 1:12,
                                     maxit = 8,
                                     measurement_model = mlpe)
{
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)
  melip.Fst <- melip.Fst[keep, keep, drop = FALSE]

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- scale_covariates(covariates)

  surface <- conductance_surface(covariates, melip.coords[keep], directions = 8)
  control <- NewtonRaphsonControl(maxit = maxit, verbose = FALSE)

  fit_linear <- suppressWarnings(
    terradish(
      melip.Fst ~ altitude + forestcover,
      data = surface,
      conductance_model = loglinear_conductance,
      measurement_model = measurement_model,
      control = control
    )
  )
  fit_spline_altitude <- suppressWarnings(
    terradish(
      melip.Fst ~ forestcover + s(altitude, df = 4),
      data = surface,
      conductance_model = smooth_loglinear_conductance,
      measurement_model = measurement_model,
      control = control
    )
  )
  fit_spline_both <- suppressWarnings(
    terradish(
      melip.Fst ~ s(altitude, df = 4) + s(forestcover, df = 4),
      data = surface,
      conductance_model = smooth_loglinear_conductance,
      measurement_model = measurement_model,
      control = control
    )
  )

  model_comparison <- data.frame(
    model = c("linear", "spline_altitude", "spline_both"),
    logLik = c(as.numeric(logLik(fit_linear)),
               as.numeric(logLik(fit_spline_altitude)),
               as.numeric(logLik(fit_spline_both))),
    AIC = c(AIC(fit_linear), AIC(fit_spline_altitude),
            AIC(fit_spline_both)),
    n_theta = c(length(coef(fit_linear)),
                length(coef(fit_spline_altitude)),
                length(coef(fit_spline_both)))
  )
  model_comparison <- model_comparison[order(model_comparison$AIC), ]

  list(
    surface = surface,
    fits = list(linear = fit_linear,
                spline_altitude = fit_spline_altitude,
                spline_both = fit_spline_both),
    model_comparison = model_comparison,
    best_model = model_comparison$model[[1]],
    best_conductance_raster = conductance(surface,
                                          list(linear = fit_linear,
                                               spline_altitude = fit_spline_altitude,
                                               spline_both = fit_spline_both)[[model_comparison$model[[1]]]]),
    best_marginal_plot = plot(
      list(linear = fit_linear,
           spline_altitude = fit_spline_altitude,
           spline_both = fit_spline_both)[[model_comparison$model[[1]]]],
      type = "marginal",
      data = surface,
      n = 40
    )
  )
}

if (sys.nframe() == 0)
{
  recovery <- run_spline_recovery_example()
  print(recovery$model_comparison)

  melip_spline <- fit_melip_spline_example()
  print(melip_spline$model_comparison)
}
