# Execute public calls before incorporating their outputs into the guides.
pkgload::load_all('.', quiet = TRUE, export_all = FALSE, compile = FALSE)
library(terra)
data(ohio_ponds, package = 'terradish')
dir.create('dev/pond_vignette_20261013', showWarnings = FALSE)
p <- ohio_ponds$ponds
groups <- factor(ohio_ponds$samples$pond, levels = p$pond)
Y <- rowsum(ohio_ponds$snp, groups, reorder = FALSE)
N <- 2 * as.numeric(table(groups))
S <- cov_from_biallelic(Y, N) / N[1]
raw <- unwrap(ohio_ponds$covariates)
fixed <- c(raw, unwrap(ohio_ponds$forest500))
surface <- conductance_surface(scale_covariates(fixed), p[c('x', 'y')], directions = 8)
raw_surface <- conductance_surface(raw[['forest']], p[c('x', 'y')], directions = 8)
control <- NewtonRaphsonControl(maxit = 200, ctol = 1e-8, ftol = 1e-8)
fit_uniform <- terradish(S ~ 1, surface, measurement_model = wishart_covariance,
                         nu = 40, optimizer = 'bfgs', control = control)
fit_forest <- terradish(S ~ forest500, surface, measurement_model = wishart_covariance,
                        nu = 40, optimizer = 'bfgs', control = control)
factory <- gaussian_smoothed_loglinear_conductance(raw_surface,
  scale_vars = 'forest', sigma_lower = 100, sigma_upper = 1500)
fit_gaussian <- terradish(S ~ forest, raw_surface, conductance_model = factory,
  measurement_model = wishart_covariance, nu = 40,
  theta = c(forest = 1, sigma.forest = 500), optimizer = 'bfgs', control = control)
fits <- list(uniform = fit_uniform, forest500 = fit_forest, gaussian = fit_gaussian)
for (nm in names(fits)) {
  x <- fits[[nm]]
  stopifnot(x$convergence$code == 0, x$fit$subproblem$convergence == 0)
  print(list(model = nm, coef = coef(x), phi = x$fit$phi,
             logLik = as.numeric(logLik(x)), AIC = AIC(x)))
}
S_ms <- cov_from_genetic_data(ohio_ponds$microsat, groups = groups,
  input = 'allele_calls', loci = ohio_ponds$microsat_loci,
  diagonal = 'gower', normalize = 'features')
fit_microsat <- terradish(S_ms ~ forest500, surface,
  measurement_model = wishart_covariance, nu = 40, optimizer = 'bfgs', control = control)
stopifnot(fit_microsat$convergence$code == 0,
          fit_microsat$fit$subproblem$convergence == 0)
print(summary(fit_forest))
print(gaussian_scale_summary(fit_gaussian))
saveRDS(list(fits = fits, fit_microsat = fit_microsat, S = S, S_ms = S_ms,
  surface = surface, raw_surface = raw_surface), 'dev/pond_vignette_20261013/smoke.rds')
print(aic_table(fits))
writeLines('complete', 'dev/pond_vignette_20261013/smoke.complete')
