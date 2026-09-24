source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/cv/setup.R")
suppressWarnings({
coords_m <- crds(pts)
S <- simulate_covariance_response(theta = c(x1 = 0.8), formula = ~x1, data = surface,
                                  tau = 1, sigma = 0.1, nu = 25, seed = 5)$covariance
D <- outer(diag(S), diag(S), "+") - 2 * S
set.seed(1); z <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("a","b","c")))
ma <- terradish(D ~ x1, surface, loglinear_conductance, mlpe_covariates(as.data.frame(z[, "a", drop = FALSE])))
mbc <- terradish(D ~ x1, surface, loglinear_conductance, mlpe_covariates(as.data.frame(z[, c("b","c")])))
cat("anova non-nested mlpe covariates:\n"); print(tryCatch(anova(ma, mbc), error = conditionMessage))
wa <- terradish(S ~ x1, surface, loglinear_conductance, wishart_covariates(as.data.frame(z[, "a", drop = FALSE]), model = "wishart_covariance"), nu = 25)
wbc <- terradish(S ~ x1, surface, loglinear_conductance, wishart_covariates(as.data.frame(z[, c("b","c")]), model = "wishart_covariance"), nu = 25)
cat("anova non-nested wishart kernels:\n"); print(tryCatch(anova(wa, wbc), error = conditionMessage))
# Gaussian-smoothed model in fixed-domain CV: baseline at theta = 0 (sigma = 0)
cov_big <- cov_r; ext(cov_big) <- c(0, 2500, 0, 2500)
pts_big <- vect(xy * 100, type = "points", crs = "local")
surf_big <- conductance_surface(cov_big, pts_big, directions = 8, saveStack = TRUE)
gfac <- gaussian_smoothed_loglinear_conductance(surf_big, scale_vars = "x1")
rf <- terradish_folds(coords_m, k = 4, method = "random", seed = 1)
cvg <- terradish_cv_folds(surf_big, list(g = D ~ x1), folds = rf, model = mlpe, conductance_model = gfac)
print(cvg$results[, c("fold","loglik","baseline_loglik","loglik_gain","error")])
})
