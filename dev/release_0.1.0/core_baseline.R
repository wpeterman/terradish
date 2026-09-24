# Run from the repository root. Every model is checkpointed immediately.
# Pass an output filename to rerun without overwriting the original baseline.
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else
  "dev/release_0.1.0/core_baseline_0.0.47.rds"
if (file.exists(out_file)) stop("Output exists; choose a new baseline filename.")
pkgload::load_all()
set.seed(20260924)
data(melip, package = "terradish")
covariates <- c(terra::unwrap(melip.altitude),
                terra::unwrap(melip.forestcover))
names(covariates) <- c("altitude", "forestcover")
coords <- terra::unwrap(melip.coords)
surface <- conductance_surface(scale_covariates(covariates), coords,
                               directions = 4, saveStack = TRUE)
control <- NewtonRaphsonControl(maxit = 200, verbose = FALSE)
records <- list(metadata = list(commit = system2("git",
  c("--no-optional-locks", "rev-parse", "HEAD"), stdout = TRUE),
  session = sessionInfo(), seed = 20260924, directions = 4,
  curvature = "exact"), models = list())
capture_fit <- function(label, formula, graph = surface,
                        factory = loglinear_conductance,
                        measurement = mlpe, nu = NULL) {
  started <- Sys.time()
  notices <- character()
  cat("START", label, format(started), "\n")
  result <- tryCatch(withCallingHandlers({
    fit <- terradish(formula, data = graph, conductance_model = factory,
      measurement_model = measurement, nu = nu,
      curvature = "exact", solver = "direct", leverage = FALSE,
      control = control)
    list(theta = coef(fit), logLik = as.numeric(logLik(fit)),
         se = sqrt(diag(-solve(fit$mle$hessian))), phi = fit$fit$phi,
         iterations = fit$cost, gradient = fit$mle$gradient,
         max_abs_gradient = max(abs(fit$mle$gradient)),
         boundary = fit$fit$boundary)
  }, warning = function(w) {
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  result$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result$warnings <- notices
  records$models[[label]] <<- result
  saveRDS(records, out_file)
  print(result)
  invisible(result)
}
capture_fit("a_mlpe", melip.Fst ~ altitude + forestcover)
capture_fit("b_leastsquares", melip.Fst ~ altitude + forestcover,
            measurement = leastsquares)
sim <- simulate_covariance_response(c(altitude = 0.3, forestcover = -0.2),
  ~ altitude + forestcover, surface, tau = 1, sigma = 0.1, nu = 100,
  seed = 20260924)
C <- sim$covariance
D <- dist_from_cov(C)
records$response <- list(C = C, D = D)
capture_fit("c_generalized_wishart", D ~ altitude + forestcover,
            measurement = generalized_wishart, nu = 100)
capture_fit("d_wishart_covariance", C ~ altitude + forestcover,
            measurement = wishart_covariance, nu = 100)
capture_fit("e_spline", melip.Fst ~ s(altitude, df = 4),
            factory = smooth_loglinear_conductance)
# A six-cell aggregate has even dimensions (26 x 30). Remove the outer
# northern/eastern strip for an odd 25 x 29 crop, retaining eligible sites.
agg <- terra::aggregate(covariates, fact = 6)
for (odd in c(TRUE, FALSE)) {
  raster <- agg
  if (odd) {
    ex <- terra::ext(raster)
    raster <- terra::crop(raster, terra::ext(ex[1], ex[2] - terra::res(raster)[1],
      ex[3], ex[4] - terra::res(raster)[2]))
  }
  site_values <- terra::extract(raster, coords, ID = FALSE)
  keep <- which(complete.cases(site_values))
  response <- melip.Fst[keep, keep]
  graph <- conductance_surface(scale_covariates(raster), coords[keep],
                               directions = 4, saveStack = TRUE)
  unit <- min(terra::res(raster))
  factory <- gaussian_smoothed_loglinear_conductance(graph,
    scale_vars = "altitude", sigma_lower = 0.5 * unit,
    sigma_upper = min(dim(raster)[1:2]) / 6 * unit)
  capture_fit(if (odd) "f_gaussian_odd" else "g_gaussian_even",
    response ~ altitude + forestcover, graph, factory)
}
site_env <- terra::extract(scale_covariates(covariates), coords, ID = FALSE)
pairwise <- pairwise_endpoint_covariates(site_env, transform = "absdiff")
capture_fit("h_mlpe_covariates", melip.Fst ~ altitude + forestcover,
            measurement = mlpe_covariates(pairwise))
kernel <- wishart_covariates(site_env["altitude"], model = "generalized_wishart")
capture_fit("i_wishart_covariates", D ~ altitude + forestcover,
            measurement = kernel, nu = 100)
stopifnot(length(records$models) == 9,
          !any(vapply(records$models, function(x) !is.null(x$error), logical(1))))
