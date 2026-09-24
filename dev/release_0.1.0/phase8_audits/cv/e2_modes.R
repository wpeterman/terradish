# Audit cv/e2: same eight noise-covariate draws, now comparing both nuisance modes.
source(file.path(root, "dev/release_0.1.0/phase8_audits/cv/setup.R"))
out <- list()
for (replicate in 1:8) {
  S <- simulate_covariance_response(c(x1 = 0.8), ~x1, surface,
    tau = 1, sigma = 0.1, nu = 25, seed = 200 + replicate)$covariance
  D <- dist_from_cov(S)
  set.seed(replicate)
  z <- rnorm(20)
  g_w <- wishart_covariates(z, model = "wishart_covariance")
  g_m <- mlpe_covariates(z)
  folds <- terradish_folds(crds(pts), k = 4, method = "random", seed = replicate)
  for (mode in c("reprofile", "fixed")) {
    for (family in c("wishart", "mlpe")) {
      response <- if (family == "wishart") S else D
      models <- if (family == "wishart") list(base = wishart_covariance, extended = g_w) else
        list(base = mlpe, extended = g_m)
      # Reprofiled comparisons are diagnostic only. The public API deliberately
      # refuses different measurement models in one reprofiled call.
      runs <- lapply(models, function(model) suppressWarnings(terradish_cv_folds(
        surface, list(m = response ~ x1), folds = folds, model = model,
        nu = if (family == "wishart") 500 else NULL, nuisance = mode)))
      base <- runs$base$results
      extended <- runs$extended$results
      stopifnot(identical(base$fold, extended$fold))
      common <- is.finite(base$loglik) & is.finite(extended$loglik)
      difference <- if (any(common)) sum(extended$loglik[common] - base$loglik[common]) else NA_real_
      row <- data.frame(replicate = replicate, family = family, nuisance = mode,
        extended_minus_base = difference, common_folds = sum(common),
        failures = sum(!is.finite(base$loglik)) + sum(!is.finite(extended$loglik)),
        baseline_failures = sum(nzchar(base$baseline_error)) + sum(nzchar(extended$baseline_error)))
      print(row)
      out[[length(out) + 1L]] <- row
    }
  }
}
comparison <- do.call(rbind, out)
print(comparison)
print(aggregate(cbind(extended_selected = comparison$extended_minus_base > 1e-6,
                      failures = comparison$failures),
  comparison[c("family", "nuisance")], sum, na.rm = TRUE))
saveRDS(comparison, file.path(root, "dev/release_0.1.0/phase8_cv_modes_comparison.rds"))
