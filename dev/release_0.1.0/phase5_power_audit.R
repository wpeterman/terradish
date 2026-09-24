.libPaths(c("dev/check/phase2-library", .libPaths()))
pkgload::load_all("dev/check/phase5-package", quiet = TRUE)
source("dev/check/phase5-package/tests/testthat/helper-robustness.R")
g <- robustness_surface()
args <- list(theta = c(x = 0.6), formula = ~x, data = g,
  sample_sizes = 6, strategies = "spacefill", nu = 50, nsim = 30,
  tau = 0.7, sigma = 0.08, seed = 905,
  control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))
matched <- do.call(covariance_response_power, args)
overconfident <- do.call(covariance_response_power, c(args, list(nu_fit = 1000)))
saveRDS(list(matched = matched, overconfident = overconfident, session = sessionInfo()),
  "dev/release_0.1.0/phase5_power_audit.rds")
print(matched$parameter_summary)
print(overconfident$parameter_summary)
print(max(abs(matched$parameter_results$estimate - overconfident$parameter_results$estimate), na.rm = TRUE))
