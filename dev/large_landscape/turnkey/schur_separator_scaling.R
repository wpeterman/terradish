#!/usr/bin/env Rscript
# schur_separator_scaling.R -- does the interface (separator) factorization
# actually become the memory ceiling of terradish_kron_reduce_tiled, and at what
# scale? For increasing raster sizes it forms the same blocks the reduction does
# and measures the Cholesky FACTOR nonzero count (memory) of (a) the largest
# per-tile interior block and (b) the assembled separator block. The separator
# factor is the candidate ceiling; this shows when it overtakes the tile factor
# and how fast it grows, so the "nested dissection at ~10M cells" premise can be
# judged on numbers.
#
# Usage: Rscript schur_separator_scaling.R --sides 256,512,724,1024,1448 \
#          --m 50 --target_tile 16000 --out schur_separator_scaling.csv
suppressMessages({library(terradish); library(terra); library(Matrix)})
# internals from the installed terradish (built from the current source)
full_laplacian <- terradish:::.terradish_full_laplacian
tile_groups    <- terradish:::.terradish_tile_groups
kron_edges     <- terradish:::.kron_edges

args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d) { i <- match(f, args); if (is.na(i)) d else args[i + 1] }
sides  <- as.integer(strsplit(ga("--sides", "256,512,724,1024,1448"), ",")[[1]])
m      <- as.integer(ga("--m", "50"))
target <- as.numeric(ga("--target_tile", "16000"))
out    <- ga("--out", "schur_separator_scaling.csv")
theta  <- c(0.6, -2.0)

factor_nnz <- function(A) {
  L <- tryCatch(Matrix::expand(Matrix::Cholesky(forceSymmetric(A), LDL = TRUE, perm = TRUE))$L,
                error = function(e) NULL)
  if (is.null(L)) NA_real_ else length(L@x)
}
mb <- function(nnz) round(nnz * 8 / 1e6, 1)   # factor values, MB

for (n in sides) {
  r  <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = n)
  rr <- terra::init(r, "row"); cc <- terra::init(r, "col")
  cov <- c((rr - n / 2) / n, (cc == (n %/% 2))); names(cov) <- c("gradient", "barrier")
  cov <- scale_covariates(cov)
  set.seed(1); surf <- conductance_surface(cov, cbind(runif(m, 1, n), runif(m, 1, n)), directions = 8)
  cond <- loglinear_conductance(~ gradient + barrier, surf$x)(theta)$conductance

  nv <- nrow(surf$x); focal <- unique(surf$demes)
  n_tiles <- max(4L, as.integer(round(nv / target)))
  groups <- tile_groups(surf, NULL, n_tiles, nv)
  tile_of <- integer(nv); for (t in seq_along(groups)) tile_of[groups[[t]]] <- t
  ep <- kron_edges(surf); cross <- tile_of[ep[, 1]] != tile_of[ep[, 2]]
  Bmask <- logical(nv); Bmask[ep[cross, 1]] <- TRUE; Bmask[ep[cross, 2]] <- TRUE
  Bmask[focal] <- TRUE; B <- which(Bmask)
  Qg <- as(full_laplacian(surf, cond), "generalMatrix")

  tile_fac <- 0; peak_int <- 0L; contribs <- list()
  for (cells in groups) {
    It <- cells[!Bmask[cells]]; if (!length(It)) next
    peak_int <- max(peak_int, length(It))
    Lii <- forceSymmetric(Qg[It, It, drop = FALSE]); tile_fac <- max(tile_fac, factor_nnz(Lii))
    LiB <- as(Qg[It, B, drop = FALSE], "CsparseMatrix"); Bt <- which(diff(LiB@p) > 0L); if (!length(Bt)) next
    Cii <- Cholesky(Lii, LDL = TRUE, perm = TRUE)
    Ct <- as.matrix(crossprod(LiB[, Bt, drop = FALSE], solve(Cii, LiB[, Bt, drop = FALSE])))
    contribs[[length(contribs) + 1L]] <- list(Bt = Bt, Ct = Ct)
  }
  S_B <- as(Qg[B, B, drop = FALSE], "generalMatrix")
  ii <- unlist(lapply(contribs, function(c) rep.int(c$Bt, length(c$Bt))), use.names = FALSE)
  jj <- unlist(lapply(contribs, function(c) rep(c$Bt, each = length(c$Bt))), use.names = FALSE)
  xx <- unlist(lapply(contribs, function(c) as.numeric(c$Ct)), use.names = FALSE)
  S_B <- as(forceSymmetric(S_B - sparseMatrix(i = ii, j = jj, x = xx, dims = dim(S_B))), "generalMatrix")
  elim_B <- setdiff(seq_along(B), match(focal, B))
  sep_fac <- factor_nnz(S_B[elim_B, elim_B, drop = FALSE])

  row <- data.frame(n = n, cells = n * n, vertices = nv, n_tiles = n_tiles,
                    peak_interior = peak_int, n_sep = length(elim_B),
                    tile_factor_nnz = tile_fac, tile_factor_MB = mb(tile_fac),
                    sep_factor_nnz = sep_fac, sep_factor_MB = mb(sep_fac),
                    sep_over_tile = round(sep_fac / tile_fac, 2))
  write.table(row, out, sep = ",", append = file.exists(out),
              col.names = !file.exists(out), row.names = FALSE)
  cat(sprintf("cells=%-8d tiles=%-3d peak_int=%-6d n_sep=%-7d | tile_fac=%sMB sep_fac=%sMB (sep/tile=%.2f)\n",
              n * n, n_tiles, peak_int, length(elim_B), mb(tile_fac), mb(sep_fac), sep_fac / tile_fac))
}
cat("done -> ", out, "\n")
