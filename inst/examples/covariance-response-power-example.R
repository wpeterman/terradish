# Worked power example for covariance-response conductance models
#
# This script demonstrates how to use covariance_response_power() to compare
# sampling designs without running forward-time population simulations.
#
# The default settings are intentionally small enough for an interactive first
# run. For study planning, increase nsim, n_designs, sample_sizes, and maxit.
#
# After installation, find this file with:
# system.file("examples", "covariance-response-power-example.R",
#             package = "terradish")

library(terradish)
library(terra)

run_covariance_response_power_example <- function(sample_sizes = c(6, 8, 10),
                                                  strategies = c("spacefill",
                                                                 "random"),
                                                  nsim = 5,
                                                  n_designs = 2,
                                                  nu = 150,
                                                  maxit = 8,
                                                  seed = 2026)
{
  data(melip, package = "terradish")
  melip.altitude <- terra::unwrap(melip.altitude)
  melip.forestcover <- terra::unwrap(melip.forestcover)
  melip.coords <- terra::unwrap(melip.coords)

  covariates <- c(melip.altitude, melip.forestcover)
  names(covariates) <- c("altitude", "forestcover")
  covariates <- scale_covariates(covariates)

  # Keep the raster stack because the Gaussian-scale model smooths raster
  # layers internally and needs access to the original grid.
  surface <- conductance_surface(
    covariates,
    melip.coords[1:12],
    directions = 8,
    saveStack = TRUE
  )

  # Define a known non-linear truth. Here forest cover has a weak negative
  # linear effect, while altitude has a hump-shaped spline effect.
  true_formula <- ~ forestcover + s(altitude, df = 4)
  true_model <- smooth_loglinear_conductance(true_formula, surface$x)
  theta_true <- attr(true_model, "default")
  theta_true["forestcover"] <- -0.15
  theta_true[grepl("^s\\(altitude\\)", names(theta_true))] <-
    c(-0.65, 1.1, -0.9, 0.45)

  # Candidate model 1: the correctly specified natural-spline model.
  #
  # Candidate model 2: a Gaussian-smoothed model that estimates the scale of
  # altitude before fitting the log-linear conductance relationship. This is not
  # the data-generating model here, but it is useful as a competing flexible
  # strategy.
  gaussian_model <- gaussian_smoothed_loglinear_conductance(
    surface,
    scale_vars = "altitude"
  )

  power <- covariance_response_power(
    theta = theta_true,
    formula = true_formula,
    data = surface,
    sample_sizes = sample_sizes,
    strategies = strategies,
    conductance_model = smooth_loglinear_conductance,
    fit_models = list(
      spline = list(
        formula = ~ forestcover + s(altitude, df = 4),
        conductance_model = smooth_loglinear_conductance
      ),
      gaussian_altitude = list(
        formula = ~ altitude + forestcover,
        conductance_model = gaussian_model,
        optimizer = "bfgs"
      ),
      linear = list(
        formula = ~ altitude + forestcover,
        conductance_model = loglinear_conductance
      )
    ),
    tau = 0.8,
    sigma = 0.12,
    nu = nu,
    nsim = nsim,
    n_designs = n_designs,
    seed = seed,
    control = NewtonRaphsonControl(maxit = maxit, verbose = FALSE),
    conductance_cor_threshold = 0.9
  )

  # Scenario-level interpretation:
  # - fit_rate should be near 1 before trusting power estimates.
  # - conductance_power is the fraction of simulations where the recovered
  #   full conductance surface correlated with truth above the chosen threshold.
  # - selected_AICc_rate is the fraction of simulations where each candidate
  #   was the best-supported model by AICc.
  scenario_summary <- power$summary[
    order(power$summary$sample_size,
          power$summary$strategy,
          power$summary$model),
    c("sample_size", "strategy", "model", "fit_rate",
      "conductance_power", "mean_conductance_cor",
      "selected_AICc_rate", "all_parameter_power")
  ]

  # Parameter-level interpretation:
  # Spline basis coefficients are useful diagnostics but are not usually the
  # main scientific target. For spline models, conductance_power and marginal
  # plots are often more interpretable than basis-coefficient power.
  parameter_summary <- power$parameter_summary[
    order(power$parameter_summary$sample_size,
          power$parameter_summary$strategy,
          power$parameter_summary$model,
          power$parameter_summary$parameter),
    c("sample_size", "strategy", "model", "parameter", "true",
      "power", "correct_sign_rate", "coverage", "bias", "rmse")
  ]

  # A simple planning table: the first sample size in each strategy/model
  # combination that reaches at least 80% conductance recovery.
  threshold <- 0.8
  planning <- scenario_summary[
    scenario_summary$fit_rate == 1 &
      scenario_summary$conductance_power >= threshold,
    c("strategy", "model", "sample_size", "conductance_power",
      "selected_AICc_rate")
  ]
  if (nrow(planning))
  {
    planning <- do.call(
      rbind,
      lapply(split(planning, list(planning$strategy, planning$model),
                   drop = TRUE),
             function(x) x[which.min(x$sample_size), , drop = FALSE])
    )
    rownames(planning) <- NULL
  }

  list(
    power = power,
    scenario_summary = scenario_summary,
    parameter_summary = parameter_summary,
    planning = planning,
    interpretation = c(
      "Start with rows where fit_rate is 1; low fit_rate means the optimizer settings need attention before interpreting power.",
      "Prefer conductance_power and mean_conductance_cor for spline recovery because basis coefficients are not directly biological parameters.",
      "Use selected_AICc_rate to ask whether the design can distinguish the flexible candidate from simpler alternatives.",
      "Increase nsim and n_designs before making final sample-size decisions."
    )
  )
}

if (sys.nframe() == 0)
{
  example <- run_covariance_response_power_example()

  cat("\nScenario summary\n")
  print(example$scenario_summary)

  cat("\nParameter summary\n")
  print(example$parameter_summary)

  cat("\nFirst sample sizes reaching 80% conductance recovery\n")
  print(example$planning)

  cat("\nInterpretation notes\n")
  writeLines(paste0("- ", example$interpretation))
}
