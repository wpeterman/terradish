.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
pkgload::load_all("dev/check/phase3-final-package")
source("tests/testthat/helper-melip.R")
output <- list()
for (n in c(8, 10)) {
  dat <- melip_fixture(seq_len(n)); S <- dat$melip.Fst
  graph <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  for (approximation in c("none", if (n == 8) "landmark" else "coarse_raster")) {
    controls <- if (approximation == "landmark") list(n_landmarks = 4, seed = 1) else list(factor = 2)
    fit <- suppressWarnings(terradish(S ~ altitude + forestcover, graph,
      measurement_model = mlpe, control = NewtonRaphsonControl(maxit = 100),
      approximation = approximation, approximation_control = controls))
    reset <- terradish:::radish_subproblem(mlpe, fit$fit$covariance, S, nu = NULL,
      control = NewtonRaphsonControl(ctol = 1e-10, ftol = 1e-10))
    print(list(n = n, approximation = approximation, theta = coef(fit),
      loglik = fit$loglik, phi = fit$fit$phi, subproblem = fit$fit$subproblem,
      reset_loglik = -reset$loglikelihood, reset_phi = reset$phi))
    output[[paste(n, approximation)]] <- list(fit = fit, reset = reset)
    saveRDS(output, "dev/release_0.1.0/warm_start_diagnosis.rds")
  }
}
