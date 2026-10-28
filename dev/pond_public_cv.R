pkgload::load_all('.', quiet = TRUE, export_all = FALSE, compile = FALSE)
data(ohio_ponds)
ponds <- ohio_ponds$ponds
g <- factor(ohio_ponds$samples$pond, levels = ponds$pond)
N <- 2 * as.numeric(table(g))
S <- cov_from_biallelic(rowsum(ohio_ponds$snp, g, reorder = FALSE), N) / N[1]
x <- c(terra::unwrap(ohio_ponds$covariates)[['forest']],
       terra::unwrap(ohio_ponds$forest500))
surface <- conductance_surface(x, ponds[c('x', 'y')], directions = 8)
factory <- gaussian_smoothed_loglinear_conductance(surface,
  scale_vars = 'forest', sigma_lower = 100, sigma_upper = 1500)
cv <- terradish_cv_folds(surface,
  formulas = list(uniform = S ~ 1, fixed500 = S ~ forest500, estimated = S ~ forest),
  folds = ponds$fold, model = wishart_covariance, nu = 40,
  conductance_model = list(uniform = loglinear_conductance,
    fixed500 = loglinear_conductance, estimated = factory),
  nuisance = 'fixed', baseline = TRUE, keep_fits = 'slim',
  control = NewtonRaphsonControl(maxit = 200, ctol = 1e-8, ftol = 1e-8),
  checkpoint = 'dev/pond_vignette_20261013/public_cv.rds')
print(cv$results)
print(cv$summary)
stopifnot(all(cv$results$error == ''),
          all(cv$results$nuisance_iterations == 0),
          all(cv$results$nuisance_convergence == 0))
write.csv(cv$results, 'dev/pond_vignette_20261013/public_cv.csv', row.names = FALSE)
writeLines('complete', 'dev/pond_vignette_20261013/public_cv.complete')
