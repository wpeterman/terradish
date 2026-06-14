#!/usr/bin/env Rscript
# schur_memory_bench.R -- demonstrate that the tiled (out-of-core) Schur/Kron
# reduction bounds peak memory on a LARGE raster, versus the single-shot
# reduction which must hold the full interior Cholesky factor.
#
# The single-shot terradish_kron_reduce() factorizes the whole interior block
# (an ~N x N sparse Laplacian); its factor fill grows superlinearly and walls on
# memory at the same scale as solver = "direct". terradish_kron_reduce_tiled()
# eliminates the interior tile by tile, so the largest dense factorization it
# ever forms is one tile's interior (~N / n_tiles), and its peak working memory
# is the sparse Laplacian plus that one tile factor -- it keeps going where the
# single-shot OOMs. Both are exact (they agree to round-off where both run).
#
# Requires an installed terradish (with terradish_kron_reduce_tiled) + terra.
# Usage:
#   Rscript schur_memory_bench.R --sides 1024,1448,2048 --m 50 \
#           --single_max_cells 1400000 --out schur_memory.csv
suppressMessages({library(terradish); library(terra); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
ga <- function(flag, d) { i <- match(flag, args); if (is.na(i)) d else args[i + 1] }
sides      <- as.integer(strsplit(ga("--sides", "1024,1448,2048"), ",")[[1]])
m          <- as.integer(ga("--m", "50"))
single_max <- as.numeric(ga("--single_max_cells", "1400000"))  # skip single-shot above this (OOM)
target     <- as.numeric(ga("--target_tile", "16000"))         # ~cells per tile
out        <- ga("--out", "schur_memory.csv")
theta      <- c(0.6, -2.0)

# tiled function is exported; fall back to ::: if running against an older build
kron_tiled <- if (exists("terradish_kron_reduce_tiled"))
  terradish_kron_reduce_tiled else terradish:::terradish_kron_reduce_tiled

build <- function(n) {
  r  <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  rr <- terra::init(r, "row"); cc <- terra::init(r, "col")
  cov <- c((rr - n / 2) / n, (cc == (n %/% 2)))
  names(cov) <- c("gradient", "barrier"); cov <- scale_covariates(cov)
  surf <- conductance_surface(cov, cbind(runif(m, 1, n), runif(m, 1, n)), directions = 8)
  list(surf = surf, cond = loglinear_conductance(~ gradient + barrier, surf$x)(theta)$conductance)
}

# realized R-heap peak (MB) and wall time around a thunk; captures errors
measure <- function(thunk) {
  invisible(gc(reset = TRUE, full = TRUE))
  t0 <- proc.time()[["elapsed"]]
  v  <- tryCatch(thunk(), error = function(e) list(.err = conditionMessage(e)))
  el <- proc.time()[["elapsed"]] - t0
  g  <- gc(full = TRUE)
  list(v = v, t = el, peakMB = round(sum(g[, 6]), 1))   # col 6 = "max used" Mb
}

# nnz of a sparse Cholesky factor (memory proxy for the single-shot interior factor)
factor_nnz <- function(A) {
  L <- tryCatch(Matrix::expand(Matrix::Cholesky(A, LDL = TRUE, perm = TRUE))$L,
                error = function(e) NULL)
  if (is.null(L)) NA_real_ else length(L@x)
}

for (n in sides) {
  N  <- as.double(n) * n
  bs <- build(n); surf <- bs$surf; cond <- bs$cond
  nv <- nrow(surf$x); interior <- setdiff(seq_len(nv), unique(surf$demes))
  nt <- max(4L, as.integer(round(nv / target)))

  # ---- single-shot (skip in the OOM regime) -------------------------------
  if (N <= single_max) {
    ss <- measure(function() terradish_kron_reduce(surf, cond))
    Q  <- terradish:::.terradish_full_laplacian(surf, cond)
    fnnz <- factor_nnz(forceSymmetric(Q[interior, interior, drop = FALSE]))
    ss_dim <- length(interior); ss_t <- round(ss$t, 2); ss_pk <- ss$peakMB
    ss_red <- if (is.null(ss$v$.err)) ss$v else NULL
    ss_status <- if (is.null(ss$v$.err)) "ok" else paste("ERR:", ss$v$.err)
  } else {
    fnnz <- NA; ss_dim <- length(interior); ss_t <- NA; ss_pk <- NA
    ss_red <- NULL; ss_status <- "skipped: > single_max_cells (OOM regime)"
  }

  # ---- tiled ---------------------------------------------------------------
  tl <- measure(function() kron_tiled(surf, cond, n_tiles = nt))
  if (is.null(tl$v$.err)) {
    tr <- tl$v; tl_status <- "ok"
    relerr <- if (!is.null(ss_red)) {
      P <- match(ss_red$boundary, tr$boundary)
      max(abs(as.matrix(ss_red$laplacian) - as.matrix(tr$laplacian[P, P]))) /
        max(abs(as.matrix(ss_red$laplacian)))
    } else NA_real_
    tl_peak_int <- tr$peak$interior; tl_peak_nnz <- tr$peak$nnz
  } else {
    tl_status <- paste("ERR:", tl$v$.err); relerr <- NA
    tl_peak_int <- NA; tl_peak_nnz <- NA
  }

  row <- data.frame(
    n = n, cells = N, vertices = nv, n_tiles = nt,
    single_dim = ss_dim, single_factor_nnz = fnnz,
    single_s = ss_t, single_peakMB = ss_pk, single_status = ss_status,
    tiled_peak_interior = tl_peak_int, tiled_peak_nnz = tl_peak_nnz,
    tiled_s = round(tl$t, 2), tiled_peakMB = tl$peakMB, tiled_status = tl_status,
    exact_relerr = signif(relerr, 3))
  write.table(row, out, sep = ",", append = file.exists(out),
              col.names = !file.exists(out), row.names = FALSE)
  cat(sprintf(paste0("n=%d cells=%.0f tiles=%d | single: dim=%d factor_nnz=%s %ss peak=%sMB [%s]",
                     " | tiled: peak_int=%s nnz=%s %ss peak=%sMB | relerr=%s\n"),
              n, N, nt, ss_dim, as.character(fnnz), as.character(ss_t), as.character(ss_pk), ss_status,
              as.character(tl_peak_int), as.character(tl_peak_nnz), round(tl$t, 2), tl$peakMB,
              as.character(signif(relerr, 3))))
}
cat("done -> ", out, "\n")
