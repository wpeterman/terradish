suppressMessages({library(terradish); library(terra); library(Matrix)})
make_landscape <- function(G = 30, n_sites = 20, seed = 1, dirs = 8)
{
  set.seed(seed)
  r <- rast(nrows = G, ncols = G, xmin = 0, xmax = G, ymin = 0, ymax = G)
  xy <- xyFromCell(r, 1:ncell(r))
  # two smooth covariates from random Fourier features
  f <- function() { k <- matrix(rnorm(20), 10); p <- runif(10, 0, 2*pi)
    z <- rowSums(sapply(1:10, function(i) cos((xy %*% k[i,]) * 0.25 + p[i]))); z }
  a <- setValues(r, f()); b <- setValues(r, f())
  cov <- c(a, b); names(cov) <- c("a", "b")
  cov <- scale_covariates(cov)
  cells <- sample(ncell(r), n_sites)
  coords <- xyFromCell(r, cells)
  surf <- conductance_surface(cov, coords, directions = dirs)
  list(cov = cov, coords = coords, surf = surf)
}
