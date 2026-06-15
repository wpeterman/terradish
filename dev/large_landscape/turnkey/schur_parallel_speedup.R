#!/usr/bin/env Rscript
# schur_parallel_speedup.R -- parallel speedup of the nested-dissection
# terradish_kron_reduce_tiled() across core counts. Recursive parallelism forks
# the two independent bisection halves (the budget halving down the recursion),
# so on a fork-capable platform the wall time should drop as cores rise WHILE
# `peak$interior` (the largest single factorization, i.e. the memory bound) stays
# flat -- the point being speed and bounded memory together.
#
# Forking is Unix-only; on Windows the nested-dissection path runs sequentially,
# so expect ~1x there. Run this on Linux/macOS (or WSL2) to see the speedup.
#
# Usage:
#   Rscript schur_parallel_speedup.R --side 1024 --m 50 --leaf 2000 \
#           --cores 1,2,4,8 --out parallel_speedup.csv
suppressMessages({library(terradish); library(terra)})
args  <- commandArgs(trailingOnly = TRUE)
ga    <- function(f, d) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
side  <- as.integer(ga("--side", "1024"))
m     <- as.integer(ga("--m", "50"))
leaf  <- as.integer(ga("--leaf", "2000"))
cores <- as.integer(strsplit(ga("--cores", "1,2,4,8"), ",")[[1]])
out   <- ga("--out", "parallel_speedup.csv")
theta <- c(0.6, -2.0)

r  <- terra::rast(nrows = side, ncols = side, xmin = 0, xmax = side, ymin = 0, ymax = side)
rr <- terra::init(r, "row"); cc <- terra::init(r, "col")
cov <- c((rr - side / 2) / side, (cc == (side %/% 2)))
names(cov) <- c("gradient", "barrier"); cov <- scale_covariates(cov)
set.seed(1)
surf <- conductance_surface(cov, cbind(runif(m, 1, side), runif(m, 1, side)), directions = 8)
cond <- loglinear_conductance(~ gradient + barrier, surf$x)(theta)$conductance
cat(sprintf("cells=%d verts=%d leaf=%d | OS=%s cores_available=%d\n",
            side * side, nrow(surf$x), leaf, .Platform$OS.type, parallel::detectCores()))

base <- NA_real_
for (k in cores) {
  t0  <- proc.time()[["elapsed"]]
  # force nested dissection: that is the path whose recursive halves parallelize
  # (the default `method = "auto"` would pick the fast single-shot here)
  res <- terradish_kron_reduce_tiled(surf, cond, method = "nested", n_tiles = leaf, cores = k)
  el  <- proc.time()[["elapsed"]] - t0
  if (is.na(base)) base <- el
  row <- data.frame(side = side, cells = side * side, leaf = leaf, cores = k,
                    seconds = round(el, 2), speedup = round(base / el, 2),
                    peak_interior = res$peak$interior)
  write.table(row, out, sep = ",", append = file.exists(out),
              col.names = !file.exists(out), row.names = FALSE)
  cat(sprintf("cores=%-2d  %6.2fs  speedup=%.2fx  peak_interior=%d\n",
              k, el, base / el, res$peak$interior))
}
cat("done -> ", out, "  (peak_interior should be constant; fork speedup is Unix-only)\n")
