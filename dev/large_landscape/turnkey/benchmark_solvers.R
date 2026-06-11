#!/usr/bin/env Rscript
# benchmark_solvers.R -- computational scaling of terradish solvers on the
# user's hardware. Builds a synthetic heterogeneous landscape at increasing
# raster sizes and times one likelihood+gradient evaluation under each solver
# (direct / amg / pcg / block_cg). No genetic simulation needed -- the response
# is generated from the model at a known theta, so this isolates solver cost.
#
# Requires an installed terradish (library(terradish) must work) + terra.
# Usage:
#   Rscript benchmark_solvers.R --sides 128,256,512,724,1024 --m 50 \
#           --solvers direct,amg,pcg --out solver_scaling.csv
#
# Notes:
#  * 'amg' is the intended large-N path (smoothed-aggregation AMG-PCG); 'direct'
#    is fastest until Cholesky fill-in exhausts memory (watch RAM at >1M cells).
#  * Increase --m to study cost vs number of focal sites at fixed N.
suppressMessages({library(terradish); library(terra)})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
sides   <- as.integer(strsplit(getarg("--sides", "128,256,512"), ",")[[1]])
m       <- as.integer(getarg("--m", "50"))
solvers <- strsplit(getarg("--solvers", "direct,amg"), ",")[[1]]
out     <- getarg("--out", "solver_scaling.csv")
theta   <- c(0.6, -2.0)
# Past direct's memory wall the Cholesky factor exhausts RAM and thrashes swap
# for a long time before aborting; cap it so direct is cleanly skipped there
# while the iterative solvers still run. Inf = no cap.
direct_max  <- as.numeric(getarg("--direct_max_cells", "Inf"))

build_surface <- function(n) {
  r  <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  rr <- terra::init(r, "row"); cc <- terra::init(r, "col")
  grad    <- (rr - n / 2) / n
  barrier <- (cc == (n %/% 2))                       # one-cell low-conductance barrier
  cov <- c(grad, barrier); names(cov) <- c("gradient", "barrier")
  cov <- scale_covariates(cov)
  pts <- cbind(runif(m, 1, n), runif(m, 1, n))
  list(surface = conductance_surface(cov, pts, directions = 8), m = m, pts = pts)
}

for (n in sides) {
  bs <- build_surface(n); surf <- bs$surface
  f  <- loglinear_conductance(~ gradient + barrier, surf$x)
  # The response only has to be a non-degenerate distance matrix for the GLS
  # starting fit: the solver work being timed comes from the internal model-E
  # solve under `solver = sv`, not from S. So synthesize S cheaply from focal
  # geometry (a noisy monotone function of Euclidean distance) and never pay for
  # a full reference solve here. Generating S from a real model solve was the
  # design flaw that hit the direct memory wall / crashed AMG before timing began.
  D  <- as.matrix(stats::dist(bs$pts)); D <- D / max(D)
  set.seed(1)
  S  <- 0.5 + 0.8 * D
  nsd   <- max(0.05 * stats::sd(S[upper.tri(S)]), 1e-3)
  noise <- matrix(stats::rnorm(bs$m * bs$m, 0, nsd), bs$m, bs$m)
  S  <- S + (noise + t(noise)) / 2; diag(S) <- 0
  for (sv in solvers) {
    if (sv == "direct" && n * n > direct_max) {
      row <- data.frame(n = n, cells = n * n, m = bs$m, solver = sv,
                        eval_s = NA, iters = NA, converged = "skipped: > direct_max_cells")
      write.table(row, out, sep = ",", append = file.exists(out),
                  col.names = !file.exists(out), row.names = FALSE)
      cat(sprintf("n=%d cells=%d solver=%-9s SKIPPED (> direct_max_cells)\n", n, n * n, sv))
      next
    }
    tm <- tryCatch({
      t0 <- proc.time()[["elapsed"]]
      fit <- terradish_algorithm(f, leastsquares, surf, S, theta, nu = 1000,
                                 gradient = TRUE, hessian = FALSE, partial = FALSE,
                                 solver = sv,
                                 solver_control = if (sv %in% c("amg","pcg","pcg_jacobi","block_cg"))
                                   list(tol = 1e-7, maxit = 5000) else list())
      list(t = proc.time()[["elapsed"]] - t0,
           iters = tryCatch({ v <- fit$solver_info$iterations
                              if (length(v)) max(v) else NA_real_ }, error = function(e) NA),
           conv  = tryCatch(all(fit$solver_info$converged), error = function(e) NA))
    }, error = function(e) list(t = NA, iters = NA, conv = paste("ERR:", conditionMessage(e))))
    row <- data.frame(n = n, cells = n * n, m = bs$m, solver = sv,
                      eval_s = round(tm$t, 3), iters = tm$iters, converged = tm$conv)
    write.table(row, out, sep = ",", append = file.exists(out),
                col.names = !file.exists(out), row.names = FALSE)
    cat(sprintf("n=%d cells=%d solver=%-9s eval=%.3fs iters=%s\n",
                n, n * n, sv, tm$t, as.character(tm$iters)))
  }
}
cat("done -> ", out, "\n")
