suppressPackageStartupMessages({library(terradish); library(terra); library(Matrix)})
make_land <- function(nr = 30, seed = 1, n_sites = 20, directions = 8) {
  set.seed(seed)
  r <- rast(nrows = nr, ncols = nr, xmin = 0, xmax = nr, ymin = 0, ymax = nr)
  # smooth random fields via focal mean of noise
  f1 <- r; values(f1) <- rnorm(ncell(r)); f1 <- focal(f1, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  f2 <- r; values(f2) <- rnorm(ncell(r)); f2 <- focal(f2, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  cov <- c(f1, f2); names(cov) <- c("x1", "x2")
  cov <- scale_covariates(cov)
  xy <- cbind(runif(n_sites, 1, nr - 1), runif(n_sites, 1, nr - 1))
  surf <- suppressWarnings(conductance_surface(cov, xy, directions = directions))
  list(surf = surf, cov = cov, xy = xy)
}
