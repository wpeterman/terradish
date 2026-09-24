suppressPackageStartupMessages({library(terradish); library(terra)})
make_land <- function(nr = 30, seed = 1, nsite = 20) {
  set.seed(seed)
  r <- rast(nrows = nr, ncols = nr, xmin = 0, xmax = nr, ymin = 0, ymax = nr)
  gx <- xFromCell(r, 1:ncell(r)); gy <- yFromCell(r, 1:ncell(r))
  # smooth-ish random covariates
  z1 <- sin(gx/5) + cos(gy/7) + rnorm(ncell(r), 0, .3)
  z2 <- (gx + gy)/nr + rnorm(ncell(r), 0, .3)
  covs <- c(setValues(r, z1), setValues(r, z2)); names(covs) <- c("x1", "x2")
  covs <- terra::focal(covs, w = 3, fun = mean, na.rm = TRUE)
  names(covs) <- c("x1", "x2")
  covs <- scale_covariates(covs)
  cells <- sample(ncell(r), nsite)
  coords <- xyFromCell(r, cells)
  surface <- conductance_surface(covs, coords, directions = 8, saveStack = TRUE)
  list(covs = covs, coords = coords, surface = surface)
}
sim_S <- function(surface, logc, nu = 500, nug = 0.2, seed = 2) {
  n <- length(surface$demes)
  cm <- terradish:::.loglinear_conductance_from_matrix(matrix(logc, ncol = 1))
  E <- as.matrix(terradish_algorithm(cm, leastsquares, surface, S = diag(n), theta = 1,
        objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)
  set.seed(seed)
  rWishart(1, df = nu, Sigma = (E + nug * diag(n)) / nu)[, , 1]
}
