# Adaptation of the supplied cv/e2.R audit. Its eight seeds are unchanged.
.libPaths(c(file.path(getwd(), "dev/check/phase2-library"), .libPaths()))
pkgload::load_all("dev/check/phase3-repaired-package")
scratch <- new.env(parent = asNamespace("terradish"))
sys.source("dev/check/phase4/R/spatial_cv.R", scratch)
mk_field <- function(seed) {
  set.seed(seed)
  r <- terra::rast(nrows = 25, ncols = 25, xmin = 0, xmax = 25, ymin = 0, ymax = 25)
  terra::values(r) <- rnorm(625)
  r <- terra::focal(r, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  (r - terra::global(r, "mean")[[1]]) / terra::global(r, "sd")[[1]]
}
rasters <- c(mk_field(1), mk_field(2)); names(rasters) <- c("x1", "x2")
set.seed(3)
coords <- terra::xyFromCell(rasters, sample(625, 20))
surface <- conductance_surface(rasters, coords, directions = 8)
output <- list()
for (replicate in 1:8) {
  S <- simulate_covariance_response(.8, ~x1, surface, tau = 1, sigma = .1,
                                     nu = 25, seed = 200 + replicate)$covariance
  D <- dist_from_cov(S)
  set.seed(replicate); z <- rnorm(20)
  folds <- terradish_folds(coords, k = 4, method = "random", seed = replicate)
  warnings <- character()
  cv <- withCallingHandlers({
    w <- scratch$terradish_cv_folds(surface, list(base = S ~ x1, env = S ~ x1), folds,
      model = list(base = wishart_covariance, env = wishart_covariates(z)),
      nu = 500, nuisance = "fixed", baseline = FALSE)
    m <- scratch$terradish_cv_folds(surface, list(base = D ~ x1, env = D ~ x1), folds,
      model = list(base = mlpe, env = mlpe_covariates(z)),
      nuisance = "fixed", baseline = FALSE)
    list(wishart = w, mlpe = m)
  }, warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  difference <- function(x) {
    scores <- setNames(x$summary$total_loglik, x$summary$model)
    difference <- unname(scores["env"] - scores["base"]); if (is.finite(difference) && abs(difference) <= 100 * .Machine$double.eps * max(1, abs(scores))) 0 else difference
  }
  row <- data.frame(replicate = replicate, wishart = difference(cv$wishart),
                     mlpe = difference(cv$mlpe),
                     failed_wishart = sum(cv$wishart$summary$failed),
                     failed_mlpe = sum(cv$mlpe$summary$failed))
  output[[replicate]] <- list(row = row, cv = cv, warnings = warnings)
  saveRDS(output, "dev/release_0.1.0/phase4_predictive_repaired.rds")
  print(row)
}
results <- do.call(rbind, lapply(output, `[[`, "row"))
print(results)
stopifnot(all(is.finite(results$wishart)), all(is.finite(results$mlpe)),
          sum(results$wishart > 0) <= 3, sum(results$mlpe > 0) <= 3)
