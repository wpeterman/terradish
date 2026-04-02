# terradish

Fast gradient-based optimization of resistance surfaces.

`terradish` is an R package for maximum likelihood estimation of isolation-by-resistance models, where conductance is a function of spatial covariates, the observed data are genetic distances, and the likelihood of the "measurement process" is cheap to compute (e.g. regression of distance matrices, or generalized Wishart). It also provides fast computation of the gradient, Hessian matrix, and derivative-based leverage/influence measures. As currently implemented it is intended for moderate-sized problems (e.g. spatial grids with less than 1mil cells, where a sparse Cholesky decomposition of the graph Laplacian is feasible). Larger problems are possible (with sufficient memory), but slow.

To put this another way: if movement across a landscape is modeled as a continuous-time Markov process, and properties of this Markov process can be empirically observed (e.g. commute time, hitting time, occupancy time; or proxies thereof), then `terradish` provides an algorithm to efficiently compute first- and second- order partial derivatives of a likelihood function with regard to entries of the infinitesimal generator of the Markov process. This is useful for fitting a parameterized Markov process -- where the entries of the infinitesimal generator depend on some low dimensional set of parameters -- to data, when the landscape is not trivially small. The algorithm implemented in `terradish` is modular, so that any twice-differentiable function of a submatrix of the inverse infinitesimal generator could conceivably be used in the likelihood.

Slides from a recent workshop can be found [here](https://github.com/nspope/radish-manuscript/raw/master/IALE_Wrkshp_Pope_Final.pdf).

![Likelihood surface for a two parameter conductance model](ms/likelihood_surface.png)

Requires [corMLPE](https://github.com/nspope/corMLPE): `devtools::install_github("nspope/corMLPE")`. Other dependencies are available through CRAN. Install `terradish` from the package source with `devtools::install()` or `remotes::install_local()`. Legacy `radish*` entry points remain available with deprecation warnings during the transition.

This is a work-in-progress and the interface is still under development. Contact at nspope at utexas dot edu.

# Worked example

```r
library(terradish)
library(terra)

data(melip)
melip.altitude <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords <- terra::unwrap(melip.coords)

# scaling spatial covariates helps avoid numeric overflow
covariates <- c(terra::scale(melip.altitude),
                terra::scale(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")

plot(covariates[["altitude"]])
points(melip.coords, pch = 19)

surface <- conductance_surface(covariates, melip.coords, directions = 8)

fit_nnls <- terradish(melip.Fst ~ forestcover + altitude, surface,
                      terradish::loglinear_conductance, terradish::leastsquares)
summary(fit_nnls)

# refit with with a different measurement model that models
# dependence among pairwise measurements (terradish::mlpe)
fit_mlpe <- terradish(melip.Fst ~ forestcover + altitude, surface,
                      terradish::loglinear_conductance, terradish::mlpe)
summary(fit_mlpe)

# visualisation:
plot(fitted(fit_mlpe, "distance"), melip.Fst, pch = 19,
     xlab = "Optimized resistance distance", ylab = "Fst")

# visualise estimated conductance surface and asymptotic confidence intervals
fitted_conductance <- conductance(surface, fit_mlpe, quantile = 0.95)

plot(fitted_conductance[["est"]], 
     main = "Fitted conductance surface\n(forestcover + altitude)")
plot(fitted_conductance[["lower95"]], 
     main = "Fitted conductance surface\n(lower 95% CI)")
plot(fitted_conductance[["upper95"]], main = 
     "Fitted conductance surface\n(upper 95% CI)")

# visualise likelihood surface across grid (takes awhile)
theta <- as.matrix(expand.grid(forestcover=seq(-1,1,length.out=21), 
                               altitude=seq(-1,1,length.out=21)))
grid <- terradish_grid(theta, melip.Fst ~ forestcover + altitude, surface,
                    terradish::loglinear_conductance, terradish::mlpe)

library(ggplot2)
ggplot(data.frame(loglik=grid$loglik, grid$theta), 
       aes(x=forestcover, y=altitude)) + 
  geom_tile(aes(fill=loglik)) + 
  geom_contour(aes(z=loglik), color="black") +
  annotate(geom = "point", colour = "red",
           x = coef(fit_mlpe)["forestcover"], 
           y = coef(fit_mlpe)["altitude"]) +
  theme_bw() +
  xlab(expression(theta[altitude])) +
  ylab(expression(theta[forestcover]))

# calculate resistance distances across grid
distances <- terradish_distance(theta, ~forestcover + altitude, 
                             surface, terradish::loglinear_conductance)

ibd <- which(theta[,1] == 0 & theta[,2] == 0)
plot(distances$distance[,,ibd], melip.Fst, pch = 19, 
     xlab = "Null resistance distance (IBD)", ylab = "Fst")

# model selection:
# fit a reduced model without "forestcover" covariate, and compare to 
# full model via a likelihood ratio test
fit_mlpe_reduced <- terradish(melip.Fst ~ altitude, surface, 
                           terradish::loglinear_conductance, terradish::mlpe)
anova(fit_mlpe, fit_mlpe_reduced)

# test for an interaction
fit_mlpe_interaction <- terradish(melip.Fst ~ forestcover * altitude, surface, 
                               terradish::loglinear_conductance, terradish::mlpe)
anova(fit_mlpe, fit_mlpe_interaction)

# test against null model of IBD
fit_mlpe_ibd <- terradish(melip.Fst ~ 1, surface,
                       terradish::loglinear_conductance, terradish::mlpe)
anova(fit_mlpe, fit_mlpe_ibd)

# concurrent IBE + IBR:
# build endpoint-difference covariates from site-level environment and
# keep the resistance surface on the conductance side
z_ibe <- pairwise_endpoint_covariates(melip.altitude, melip.coords,
                                      transform = "absdiff", scale = TRUE)
g_ibe <- mlpe_covariates(z_ibe)

fit_ibe_only <- terradish(melip.Fst ~ 1, surface,
                       terradish::loglinear_conductance, g_ibe)
fit_ibr_only <- fit_mlpe
fit_joint <- terradish(melip.Fst ~ forestcover + altitude, surface,
                    terradish::loglinear_conductance, g_ibe)

# compare IBD, IBE-only, IBR-only, and joint models
cv_model_selection(list(
  list(train_mod = fit_mlpe_ibd, cv_loglik = fit_mlpe_ibd$loglik, full_mod = fit_mlpe_ibd),
  list(train_mod = fit_ibe_only, cv_loglik = fit_ibe_only$loglik, full_mod = fit_ibe_only),
  list(train_mod = fit_ibr_only, cv_loglik = fit_ibr_only$loglik, full_mod = fit_ibr_only),
  list(train_mod = fit_joint, cv_loglik = fit_joint$loglik, full_mod = fit_joint)
), aic = TRUE)

# if you already have fitted terradish models and only want an
# information-criterion table, use aic_table()
aic_table(
  list(fit_mlpe_ibd, fit_ibe_only, fit_ibr_only, fit_joint),
  mod_names = c("IBD", "IBE", "IBR", "IBE + IBR"),
  AICc = TRUE
)

# for predictive assessment on repeated random train/test splits,
# use terradish_cv_replicates()
cv_reps <- terradish_cv_replicates(
  melip.coords,
  covariates,
  melip.Fst ~ forestcover + altitude,
  model = terradish::mlpe,
  seeds = c(1, 2, 3),
  fit_full = FALSE
)
summary(cv_reps)

# categorical covariates:
# categorical raster layers should be factor-valued, see ?terra::as.factor
# the names of levels are taken from the VALUE column when it exists
forestcover_class <- cut(terra::values(melip.forestcover)[,1], breaks = c(0, 1/6, 1/3, 1))
melip.forestcover_cat <- terra::setValues(melip.forestcover, as.numeric(forestcover_class))
melip.forestcover_cat <- terra::as.factor(melip.forestcover_cat)

RAT <- levels(melip.forestcover_cat)[[1]]
RAT$VALUE <- levels(forestcover_class) #explicitly defines names in the RAT
levels(melip.forestcover_cat) <- RAT

covariates_cat <- c(melip.forestcover_cat, melip.altitude)
names(covariates_cat) <- c("forestcover", "altitude")

surface_cat <- conductance_surface(covariates_cat, melip.coords, directions = 8)

fit_mlpe_cat <- terradish(melip.Fst ~ forestcover + altitude, surface_cat, 
                       terradish::loglinear_conductance, terradish::mlpe)

# contrast coding is the default for R, and for this conductance model
# the (non-identifiable) intercept is omitted (e.g. only relative
# differences in conductance among levels are identifiable from the data)
summary(fit_mlpe_cat) 

# example of lower level interface:
# compute negative loglikelihood, gradient, Hessian for a given choice of
# of the conductance parameters theta, using a different measurement model
# (terradish::generalized_wishart)
terradish_algorithm(terradish::loglinear_conductance(~forestcover + altitude, surface$x), 
                 terradish::generalized_wishart, surface, 
                 pmax(melip.Fst, 0), nu = 1000, theta = c(-0.3, 0.3), 
                 gradient = TRUE, hessian = TRUE)$hessian
# numerical verification (not run)
#numDeriv::hessian(function(x)
#     terradish_algorithm(terradish::loglinear_conductance(~forestcover + altitude, surface$x), 
#                      terradish::generalized_wishart, surface, 
#                      pmax(melip.Fst, 0), nu = 1000, theta = x)$objective,
#                  c(-0.3, 0.3))
```

# Helper workflows

In addition to the core fitting functions, `terradish` includes a few helpers
that are useful once you start comparing models or organizing larger analyses:

- `aic_table()` ranks fitted `terradish` models by AIC, AICc, or BIC when the
  models were fit to the same focal set and landscape surface.
- `cv_model_selection()` is helpful when you already have cross-validation
  results and want to compare held-out loglikelihood together with an
  information-criterion table from the full-data fits.
- `terradish_cv_replicates()` repeats the same cross-validation workflow across
  multiple random splits, which is often more informative than relying on a
  single train/test partition.
- `terradish_results()` inspects a saved terradish-style results directory, and
  `terradish_parameters()` extracts a compact coefficient table from a saved
  fitted model.

The vignette `vignette("ibe-ibr-workflow", package = "terradish")` gives a
more guided walkthrough of these helper functions with additional context and
annotated examples.
 
# RStan hooks
In progress
