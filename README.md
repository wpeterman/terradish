---
editor_options:
  markdown:
    wrap: 72
---

# terradish

Fast gradient-based optimization of resistance surfaces.

`terradish` is an R package for maximum likelihood estimation of
isolation-by-resistance models, where conductance is a function of
spatial covariates, the observed data are genetic distances or
covariance summaries, and the measurement model remains cheap enough to
profile efficiently. It provides fast derivatives for conductance-model
parameters, supports several likelihood layers, and is designed for
moderate-sized raster problems where sparse Laplacian factorization is
still feasible.

The package is now `terra`-native and includes both the core
`terradish()` fitting workflow and helper tools for model comparison,
repeated cross-validation, saved-results extraction, covariance-based
simulation, and Gaussian scale-of-effect fitting.

To put this another way: if movement across a landscape is modeled as a
continuous-time Markov process, and properties of this Markov process
can be empirically observed through pairwise genetic distances or
related covariance summaries, then `terradish` provides an efficient way
to fit low-dimensional conductance models to those data.

Useful entry points include:

-   `conductance_surface()` to build a graph from `terra` rasters and
    focal locations
-   `terradish()` to fit conductance models with `leastsquares`, `mlpe`,
    `generalized_wishart`, or `wishart_covariance`
-   `plot(..., type = "fit" | "marginal" | "marginal_response")` for
    fitted-model visualization
-   `gaussian_smoothed_loglinear_conductance()` for joint optimization
    of conductance coefficients and Gaussian scale-of-effect parameters
-   `aic_table()`, `cv_model_selection()`, and
    `terradish_cv_replicates()` for comparing fitted models
-   `simulate_covariance_response()` for generating covariance responses
    under the `wishart_covariance` measurement model


Most dependencies are available through CRAN. Install `terradish` from
package source with `devtools::install()` or `remotes::install_local()`.
Legacy `radish*` entry points remain available with deprecation warnings
during the transition.

This is still an active development package. Contact Bill Peterman
(Peterman.73\@osu.edu) or submit an issue on GitHub if you encounter a
problem or have a feature request.

# Worked Example

``` r
library(terradish)
library(terra)

data(melip)
melip.altitude <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords <- terra::unwrap(melip.coords)

# scaling spatial covariates helps avoid numeric overflow
covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

plot(covariates[["altitude"]])
points(melip.coords, pch = 19)

surface <- conductance_surface(covariates, melip.coords, directions = 8)

fit_nnls <- terradish(
  melip.Fst ~ forestcover + altitude,
  surface,
  terradish::loglinear_conductance,
  terradish::leastsquares
)
summary(fit_nnls)

# refit with a measurement model that accounts for dependence among
# pairwise measurements
fit_mlpe <- terradish(
  melip.Fst ~ forestcover + altitude,
  surface,
  terradish::loglinear_conductance,
  terradish::mlpe
)
summary(fit_mlpe)

# fitted-value plot
plot(fit_mlpe, type = "fit")

# fitted conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

plot(
  fitted_conductance[["est"]],
  main = "Fitted conductance surface\n(forestcover + altitude)"
)
plot(
  fitted_conductance[["lower95"]],
  main = "Fitted conductance surface\n(lower 95% CI)"
)
plot(
  fitted_conductance[["upper95"]],
  main = "Fitted conductance surface\n(upper 95% CI)"
)

# marginal effects on the conductance and response scales
plot(fit_mlpe, type = "marginal", data = surface)
plot(fit_mlpe, type = "marginal_response", data = surface)

# likelihood surface across a parameter grid
theta <- as.matrix(
  expand.grid(
    forestcover = seq(-1, 1, length.out = 21),
    altitude = seq(-1, 1, length.out = 21)
  )
)
grid <- terradish_grid(
  theta,
  melip.Fst ~ forestcover + altitude,
  surface,
  terradish::loglinear_conductance,
  terradish::mlpe
)

library(ggplot2)
ggplot(data.frame(loglik = grid$loglik, grid$theta),
       aes(x = forestcover, y = altitude)) +
  geom_tile(aes(fill = loglik)) +
  geom_contour(aes(z = loglik), color = "black") +
  annotate(
    geom = "point",
    colour = "red",
    x = coef(fit_mlpe)["forestcover"],
    y = coef(fit_mlpe)["altitude"]
  ) +
  theme_bw() +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))
```

# Gaussian Scale Optimization

`terradish` can now optimize Gaussian scale-of-effect parameters inside
the conductance model rather than treating raster smoothing as an outer
pre-processing step.

``` r
surface_raw <- conductance_surface(
  melip.forestcover,
  melip.coords,
  directions = 8,
  saveStack = TRUE
)

fit_gaussian <- terradish(
  melip.Fst ~ forestcover,
  data = surface_raw,
  conductance_model = gaussian_smoothed_loglinear_conductance(surface_raw),
  measurement_model = terradish::mlpe,
  optimizer = "auto",
  leverage = FALSE
)

gaussian_scale_summary(fit_gaussian)
plot(fit_gaussian, type = "sigma")
```

The vignette
`vignette("gaussian-scale-optimization", package = "terradish")` gives a
staged workflow for fitting one-raster scale-aware models, inspecting
`sigma`, comparing them to fixed-raster fits, and only then adding
additional landscape variables.

# Helper Workflows

In addition to the core fitting functions, `terradish` includes helper
functions that are useful once you start comparing models or organizing
larger analyses:

-   `aic_table()` ranks fitted `terradish` models by AIC, AICc, or BIC
    when the models were fit to the same focal set and landscape surface
-   `cv_model_selection()` compares held-out log-likelihood together
    with information-criterion summaries from the corresponding
    full-data fits
-   `terradish_cv_replicates()` repeats the same train/test workflow
    across multiple random splits
-   `terradish_results()` inspects a saved terradish-style results
    directory, and `terradish_parameters()` extracts a compact parameter
    table from saved models
-   `simulate_covariance_response()` generates covariance responses from
    a known conductance surface under the `wishart_covariance`
    measurement model

The vignette `vignette("ibe-ibr-workflow", package = "terradish")`
provides a more guided walkthrough of these helper functions with
additional context and annotated examples.
