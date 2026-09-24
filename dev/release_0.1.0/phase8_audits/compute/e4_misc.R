source("C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/terradish/dev/release_0.1.0/phase8_audits/compute/common.R")
cat("==== E4b crop_to_focal_buffer boundary effect ====\n")
set.seed(5)
G <- 40
r <- rast(nrows = G, ncols = G, xmin = 0, xmax = G, ymin = 0, ymax = G)
xy <- xyFromCell(r, 1:ncell(r))
k <- matrix(rnorm(20), 10); p <- runif(10, 0, 2*pi)
z <- rowSums(sapply(1:10, function(i) cos((xy %*% k[i,]) * 0.25 + p[i])))
cv <- setValues(r, z); names(cv) <- "a"; cv <- scale_covariates(cv)
cells <- cellFromRowCol(r, rep(c(14, 20, 26), each = 3), rep(c(14, 20, 26), 3))
pts <- xyFromCell(r, cells)
Rfull <- NULL
for (buf in c(NA, 0, 2, 5, 10)) {
  ss <- if (is.na(buf)) conductance_surface(cv, pts, directions = 8) else conductance_surface(cv, pts, directions = 8, crop_buffer = buf)
  for (th in c(0, 1.5)) {
    m <- loglinear_conductance(~ a, ss$x)
    rr <- terradish_distance(matrix(th, 1, dimnames = list(NULL, "a")), ~ a, ss, loglinear_conductance)$distance[, , 1]
    if (is.na(buf)) { if (th == 0) R0 <- rr else R1 <- rr; next }
    ref <- if (th == 0) R0 else R1
    off <- upper.tri(rr)
    cat(sprintf("buffer=%2d cells (vertices %4d) theta=%.1f: mean rel increase in R = %.3f, max = %.3f, cor with full = %.4f\n",
                buf, nrow(ss$x), th, mean(rr[off]/ref[off] - 1), max(rr[off]/ref[off] - 1), cor(rr[off], ref[off])))
  }
}
# coefficient shift: simulate on full, fit on cropped
sfull <- conductance_surface(cv, pts, directions = 8)
simc <- simulate_covariance_response(theta = c(a = 1), formula = ~ a, data = sfull, tau = 1, sigma = 0.02, nu = 500, seed = 9)
Sc <- simc$covariance
for (buf in c(NA, 2, 5, 10)) {
  ss <- if (is.na(buf)) sfull else conductance_surface(cv, pts, directions = 8, crop_buffer = buf)
  fit <- suppressWarnings(terradish(Sc ~ a, ss, measurement_model = wishart_covariance, nu = 500))
  cat(sprintf("fit on %s: theta_a = %.3f  loglik = %.2f\n", ifelse(is.na(buf), "full", paste("buffer", buf)), coef(fit), fit$loglik))
}
