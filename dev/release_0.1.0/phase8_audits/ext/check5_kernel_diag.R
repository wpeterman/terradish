# Does the new default absdiff kernel respond to diagonal heterogeneity aligned with z^2
# (no environment-aligned off-diagonal covariance in the truth)?
suppressMessages({library(terradish); library(terra)})
data(melip)
alt <- unwrap(melip.altitude); fc <- unwrap(melip.forestcover); co <- unwrap(melip.coords)
cov <- aggregate(c(alt, fc), fact = 5, na.rm = TRUE)
names(cov) <- c("altitude", "forestcover")
cov <- scale_covariates(cov)
keep <- 1:20
surface <- conductance_surface(cov, co[keep], directions = 8)
ctrl <- NewtonRaphsonControl(maxit = 60, verbose = FALSE)
n <- 20
f <- loglinear_conductance(~ altitude, surface$x)
E <- as.matrix(terradish_algorithm(f, leastsquares, surface, S = diag(n), theta = 0.5,
                                   objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
set.seed(3)
z <- rnorm(n); zc <- z - mean(z)
nug <- exp(log(0.2) + 0.5 * zc^2)     # diagonal-only heterogeneity, larger at env extremes
Sigma <- E + diag(nug)
nu <- 200
out <- NULL
for (s in 1:4) {
  set.seed(100 + s)
  S <- rWishart(1, nu, Sigma)[, , 1] / nu
  D <- outer(diag(S), rep(1, n)) + outer(rep(1, n), diag(S)) - 2 * S
  for (mdl in c("wishart_covariance", "generalized_wishart")) {
    base <- if (mdl == "wishart_covariance") wishart_covariance else generalized_wishart
    Y <- if (mdl == "wishart_covariance") S else D
    f0 <- suppressWarnings(terradish(Y ~ altitude, surface, loglinear_conductance, base, nu = nu, control = ctrl, leverage = FALSE))
    f1 <- suppressWarnings(terradish(Y ~ altitude, surface, loglinear_conductance,
                                     wishart_covariates(data.frame(env = z), model = mdl, transform = "absdiff"), nu = nu, control = ctrl, leverage = FALSE))
    out <- rbind(out, data.frame(seed = s, model = mdl, lambda = signif(f1$fit$phi["lambda_absdiff_env", 1], 3),
                                 dAIC = signif(f1$aic - f0$aic, 4)))
  }
}
print(out)
