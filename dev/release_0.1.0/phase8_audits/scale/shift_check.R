suppressMessages({library(terra); library(terradish)})
prep_fun <- terradish:::.gaussian_scale_prepare_layer
lv <- terradish:::.gaussian_scale_layer_values
brute <- function(m, sig) {  # normalized Gaussian smoothing, zero padding, unit cells
  nr <- nrow(m); nc <- ncol(m); out <- m
  for (i in 1:nr) for (j in 1:nc) {
    w <- outer((1:nr - i)^2, (1:nc - j)^2, `+`); w <- exp(-w / (2 * sig^2))
    out[i, j] <- sum(w * m) / sum(w)
  }
  out
}
set.seed(3)
for (dims in list(c(30, 30), c(31, 31), c(30, 31), c(31, 30))) {
  nr <- dims[1]; nc <- dims[2]
  r <- rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc, ymin = 0, ymax = nr); values(r) <- rnorm(nr * nc)
  rc <- as.matrix(expand.grid(row = 1:nr, col = 1:nc))
  v <- lv(prep_fun(r, rc), 2, standardize = FALSE)$value
  m <- matrix(NA, nr, nc); m[rc] <- v
  b <- brute(as.matrix(r, wide = TRUE), 2)
  inner_r <- 3:(nr - 3); inner_c <- 3:(nc - 3)
  d0 <- max(abs(m[inner_r, inner_c] - b[inner_r, inner_c]))
  # test shifted alignment: package[i,j] == brute[i+dr, j+dc]
  best <- NULL
  for (dr in -1:1) for (dc in -1:1) {
    d <- max(abs(m[inner_r, inner_c] - b[inner_r + dr, inner_c + dc]))
    if (is.null(best) || d < best[3]) best <- c(dr, dc, d)
  }
  cat(sprintf("%dx%d: max|pkg - brute| (interior) = %.2e; best alignment pkg[i,j]=brute[i%+d,j%+d] with err %.2e\n",
              nr, nc, d0, best[1], best[2], best[3]))
}
# Also check that the kernel window truncation at +-half extent affects values for large sigma
