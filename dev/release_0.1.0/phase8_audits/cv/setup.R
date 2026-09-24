suppressMessages({library(terradish); library(terra)})
set.seed(11)
n <- 25
mk_field <- function(seed) {
  set.seed(seed)
  r <- rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  values(r) <- rnorm(n * n)
  r <- focal(r, w = matrix(1, 5, 5), fun = mean, na.rm = TRUE)
  r <- (r - global(r, "mean")[[1]]) / global(r, "sd")[[1]]
  r
}
x1 <- mk_field(1); x2 <- mk_field(2)
cov_r <- c(x1, x2); names(cov_r) <- c("x1", "x2")
set.seed(3)
cells <- sample(ncell(x1), 20)
xy <- xyFromCell(x1, cells)
pts <- vect(xy, type = "points", crs = "local")
surface <- conductance_surface(cov_r, pts, directions = 8)
