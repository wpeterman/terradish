source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/core/helper.R")
set.seed(2)
nr <- 30
r <- rast(nrows = nr, ncols = nr, xmin = 0, xmax = nr, ymin = 0, ymax = nr)
f1 <- r; values(f1) <- rnorm(ncell(r)); f1 <- focal(f1, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
f2 <- r; values(f2) <- rnorm(ncell(r)); f2 <- focal(f2, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
cov <- c(f1, f2); names(cov) <- c("x1", "x2"); cov <- scale_covariates(cov)
xy <- cbind(runif(20, 1, 12), runif(20, 1, 12))   # sites clustered in one corner
surf <- suppressWarnings(conductance_surface(cov, xy, directions = 8))
theta_true <- c(x1 = 0.6, x2 = -0.4); tau_true <- 2; sigma_true <- 0.5; nu <- 300
n <- length(surf$demes); H <- diag(n) - 1/n
res <- list()
for (rep in 1:5) {
  sim <- simulate_covariance_response(theta_true, ~ x1 + x2, surf, tau = tau_true,
                                      sigma = sigma_true, nu = nu, seed = 100 + rep)
  C <- sim$covariance
  if (rep == 1) { E <- sim$E; cat("1'E1/n^2 =", sum(E)/n^2, " mean diag(E) =", mean(diag(E)), "\n") }
  Cc <- H %*% C %*% H
  fit <- function(S, g) suppressWarnings(terradish(S ~ x1 + x2, surf, loglinear_conductance, g, nu = nu))
  f_wc <- fit(C, wishart_covariance); f_wcc <- fit(Cc, wishart_covariance); f_gw <- fit(dist_from_cov(C), generalized_wishart)
  row <- function(f, lab) c(model = lab, rep = rep, round(c(coef(f), tau = f$fit$phi[1], sigma = exp(f$fit$phi[2]),
        se_x1 = tryCatch(sqrt(diag(solve(f$fit$hessian)))[1], error=function(e) NA), bnd = f$fit$boundary), 4))
  res[[length(res)+1]] <- row(f_wc, "wc_uncentered"); res[[length(res)+1]] <- row(f_wcc, "wc_centered"); res[[length(res)+1]] <- row(f_gw, "gw_dist")
}
tabn <- as.data.frame(do.call(rbind, res)); for (j in 3:ncol(tabn)) tabn[[j]] <- as.numeric(tabn[[j]])
print(tabn)
print(aggregate(tabn[, 3:ncol(tabn)], list(model = tabn$model), mean))
