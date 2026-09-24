# Audit cv/e3 converted from the removed split API to fixed-domain folds.
source(file.path(root, "dev/release_0.1.0/phase8_audits/cv/setup.R"))
cov_big <- cov_r
ext(cov_big) <- c(0, 2500, 0, 2500)
pts_big <- vect(xy * 100, type = "points", crs = "local")
surf_big <- conductance_surface(cov_big, pts_big, directions = 8, saveStack = TRUE)
gfac <- gaussian_smoothed_loglinear_conductance(surf_big, scale_vars = "x1")
S <- simulate_covariance_response(c(x1 = 0.8), ~x1, surf_big,
  tau = 1, sigma = 0.1, nu = 25, seed = 5)$covariance
D <- dist_from_cov(S)
folds <- terradish_folds(xy, k = 4, method = "random", seed = 1)
cv <- terradish_cv_folds(surf_big, list(g = D ~ x1), folds = folds,
  model = mlpe, conductance_model = gfac, nuisance = "reprofile", keep_fits = "full")
print(cv$results[, c("fold", "loglik", "baseline_loglik", "error", "baseline_error")])
differences <- numeric()
for (i in which(is.finite(cv$results$loglik))) {
  row <- cv$results[i, ]
  fit <- cv$fits[[paste("g", 1, row$fold, sep = "::")]]
  test <- which(folds == as.integer(row$fold))
  graph <- terradish:::.terradish_subset_graph_focal(surf_big, test)
  f <- fit$submodels$f_factory(~x1, graph$x)
  score <- function(theta) -terradish_algorithm(f, mlpe, graph,
    D[test, test, drop = FALSE], theta = theta, gradient = FALSE,
    hessian = FALSE, partial = FALSE)$objective
  correct <- score(fit$mle$theta_internal)
  wrong <- score(coef(fit))
  differences <- c(differences, correct - row$loglik)
  cat("Fold", row$fold, "CV", row$loglik, "internal", correct,
      "incorrect external-scale input", wrong, "\n")
}
stopifnot(length(differences) > 0, max(abs(differences)) < 1e-6)
