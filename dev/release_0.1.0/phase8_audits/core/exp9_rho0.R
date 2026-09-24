source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/core/helper.R")
L <- make_land(25, seed = 4, n_sites = 15)
surf <- L$surf
R <- terradish_distance(theta = matrix(c(0.7, -0.3), 1), formula = ~x1 + x2, data = surf, conductance_model = loglinear_conductance)$distance[, , 1]
set.seed(3); n <- nrow(R)
S <- 0.01 + 0.02 * R + { e <- matrix(rnorm(n * n, 0, 0.002), n); e[upper.tri(e)] <- t(e)[upper.tri(e)]; e }; diag(S) <- 0
fm <- withCallingHandlers(terradish(S ~ x1 + x2, surf, loglinear_conductance, mlpe), warning = function(w) {cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")})
fl <- terradish(S ~ x1 + x2, surf, loglinear_conductance, leastsquares)
cat("mlpe phi:", signif(fm$fit$phi, 4), " rho_actual =", signif(plogis(fm$fit$phi[4])/2, 3), "\n")
cat("mlpe theta/SE:", signif(coef(fm), 4), signif(sqrt(diag(solve(fm$fit$hessian))), 4), "\n")
cat("LS   theta/SE:", signif(coef(fl), 4), signif(sqrt(diag(solve(fl$fit$hessian))), 4), "\n")
print(summary(fm)$phi_table)
cat("subproblem:", unlist(fm$fit$subproblem), "\n")
cat("Outer convergence record:\n"); print(fm$convergence)
stopifnot(max(abs(coef(fm) - coef(fl))) < 1e-5,
          max(abs(sqrt(diag(vcov(fm))) - sqrt(diag(vcov(fl))))) < 1e-5,
          summary(fm)$phi_table["rho", "Estimate"] >= 0,
          summary(fm)$phi_table["rho", "Estimate"] < 1e-6)
