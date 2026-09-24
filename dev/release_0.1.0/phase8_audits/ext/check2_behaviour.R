suppressMessages({library(terradish); library(terra)})
data(melip)
alt <- unwrap(melip.altitude); fc <- unwrap(melip.forestcover); co <- unwrap(melip.coords)
cov <- aggregate(c(alt, fc), fact = 5, na.rm = TRUE)
names(cov) <- c("altitude", "forestcover")
cov <- scale_covariates(cov)
cat("raster dim:", dim(cov)[1:2], "\n")
keep <- 1:20
co20 <- co[keep]
surface <- conductance_surface(cov, co20, directions = 8)
ctrl <- NewtonRaphsonControl(maxit = 50, verbose = FALSE)

cat("\n==== (j) `scale changes absdiff kernel magnitude ====\n")
site_alt <- terra::extract(alt, co20)[, 2]
k1 <- attr(wishart_covariates(site_alt, scale = FALSE), "kernel_covariates")
k2 <- attr(wishart_covariates(site_alt, scale = TRUE), "kernel_covariates")
cat("identical kernels scale=FALSE vs TRUE:", isTRUE(all.equal(unclass(k1)[, , 1], unclass(k2)[, , 1], check.attributes = FALSE)), "\n")

cat("\n==== (k) absdiff kernel uses centered pairwise differences ====\n")
zz <- c(-2, -2, 0, 0, 2, 2)
K <- attr(wishart_covariates(zz), "kernel_covariates")[, , 1]
print(round(K, 2))
cat("sites 3,4 identical env (at mean): K34 =", K[3, 4], "; sites 5,6 identical env (extreme): K56 =", K[5, 6],
    "; sites 1,5 opposite: K15 =", K[1, 5], "\n")
cat("implied squared-distance contribution K_ii+K_jj-2K_ij equals unscaled abs difference:",
    isTRUE(all.equal(outer(diag(K), rep(1, 6)) + outer(rep(1, 6), diag(K)) - 2 * K,
                     abs(outer(zz, zz, "-")), check.attributes = FALSE)), "\n")

cat("\n==== (f) boundary LRT for lambda: anova now uses chi-bar-square ====\n")
nu <- 200
res <- list()
for (s in 1:4) {
  sim <- simulate_covariance_response(theta = c(altitude = 0.4, forestcover = -0.3),
                                      formula = ~ altitude + forestcover, data = surface,
                                      tau = 1, sigma = 0.2, nu = nu, seed = s)
  Sc <- sim$covariance
  env <- rnorm(length(keep))   # pure-noise environmental covariate
  f0 <- suppressWarnings(terradish(Sc ~ altitude + forestcover, surface, loglinear_conductance,
                                   wishart_covariance, nu = nu, control = ctrl))
  f1 <- suppressWarnings(terradish(Sc ~ altitude + forestcover, surface, loglinear_conductance,
                                   wishart_covariates(env, model = "wishart_covariance"), nu = nu, control = ctrl))
  a <- tryCatch(anova(f1, f0), error = function(e) conditionMessage(e))
  lam <- f1$fit$phi["lambda_absdiff_var1", 1]
  pt <- summary(f1)$phi_table["lambda_absdiff_var1", ]
  res[[s]] <- data.frame(seed = s, lambda = signif(lam, 3), LR = if (is.character(a)) NA else signif(a[2, "ChiSq"], 3),
                         p_anova = if (is.character(a)) a else signif(a[2, "Pr(>Chi)"], 3),
                         p_chibar = if (is.character(a)) NA else signif(0.5 * pchisq(a[2, "ChiSq"], 1, lower.tail = FALSE), 3),
                         wald_lo = signif(pt[3], 3), wald_hi = signif(pt[4], 3),
                         dAIC = signif(f1$aic - f0$aic, 3), boundary_flag = f1$fit$boundary)
}
print(do.call(rbind, res))

cat("\n==== (g) anova must reject non-nested measurement models ====\n")
Fst <- melip.Fst[keep, keep]
site <- data.frame(alt = site_alt, fc = terra::extract(fc, co20)[, 2])
xy <- terra::geom(co20)[, c("x", "y")]
gA <- mlpe_covariates(site["alt"], transform = "absdiff", scale = TRUE)
gB <- mlpe_covariates(cbind(site["fc"], x = xy[, 1]), transform = "absdiff", scale = TRUE)
fA <- suppressWarnings(terradish(Fst ~ altitude, surface, loglinear_conductance, gA, control = ctrl))
fB <- suppressWarnings(terradish(Fst ~ altitude, surface, loglinear_conductance, gB, control = ctrl))
print(tryCatch(anova(fB, fA), error = function(e) paste("ERROR:", conditionMessage(e))))


cat("Corrected README argument order:"); print(dim(pairwise_endpoint_covariates(cov, co20)))
