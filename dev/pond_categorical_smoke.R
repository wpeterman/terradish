pkgload::load_all('.', quiet = TRUE, export_all = FALSE, compile = FALSE)
data(ohio_ponds)
library(terra)
cover <- unwrap(readRDS('inst/extdata/ohio-ponds/landcover100.rds'))
slope <- scale_covariates(unwrap(ohio_ponds$covariates)[['slope']])
surface <- conductance_surface(c(cover, slope), ohio_ponds$ponds[c('x', 'y')])
surface$x$landcover <- relevel(surface$x$landcover, ref = 'forest')
contrasts(surface$x$landcover) <- contr.treatment(levels(surface$x$landcover), base = 1)
g <- factor(ohio_ponds$samples$pond, levels = ohio_ponds$ponds$pond)
N <- 2 * as.numeric(table(g))
S <- cov_from_biallelic(rowsum(ohio_ponds$snp, g, reorder = FALSE), N) / N[1]
fit <- terradish(S ~ landcover, surface, measurement_model = wishart_covariance,
  nu = 40, optimizer = 'bfgs',
  control = NewtonRaphsonControl(maxit = 200, ctol = 1e-8, ftol = 1e-8))
print(summary(fit))
stopifnot(fit$convergence$code == 0, fit$fit$subproblem$convergence == 0)
print(exp(coef(fit)))
cv <- terradish_cv_folds(surface,
  list(uniform = S ~ 1, categorical = S ~ landcover, mixed = S ~ landcover + slope),
  folds = ohio_ponds$ponds$fold, model = wishart_covariance, nu = 40,
  nuisance = 'fixed', optimizer = 'bfgs', keep_fits = 'none',
  control = NewtonRaphsonControl(maxit = 200, ctol = 1e-8, ftol = 1e-8),
  checkpoint = 'dev/pond_vignette_20261013/categorical_cv.rds')
print(cv$summary)
stopifnot(all(cv$results$error == ''))
write.csv(cv$results, 'dev/pond_vignette_20261013/categorical_cv.csv', row.names = FALSE)
writeLines('complete', 'dev/pond_vignette_20261013/categorical.complete')
