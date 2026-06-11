# Computational scaling benchmark: one raster size per invocation.
# Builds an 8-neighbour grid Laplacian with heterogeneous conductance, a block
# RHS of m focal sites, and times the solve via direct Cholesky, IC-PCG (Eigen
# IncompleteCholesky CG, compiled), and block-CG. Records solve time, iterations,
# direct-factor fill-in (a memory proxy), and accuracy vs the direct solve.
# Appends one row to a CSV so a sweep over N can run across calls.
#
# Run: TERRADISH_SO=/tmp/tdsrc/src/terradish.so Rscript bench_scaling.R <n_side> <m> <out.csv>
suppressMessages({library(Matrix); library(Rcpp)})
PKG <- Sys.getenv("TERRADISH_PKG", "/sessions/funny-vibrant-bohr/mnt/terradish")
SO  <- Sys.getenv("TERRADISH_SO",  "/tmp/tdsrc/src/terradish.so")
dll <- dyn.load(SO)
for (s in c("assemble_reduced_laplacian","pcg_reduced_laplacian_ic",
            "block_cg_reduced_laplacian"))
  assign(paste0("_terradish_", s), getNativeSymbolInfo(paste0("_terradish_", s), dll), envir = .GlobalEnv)
source(file.path(PKG, "R", "RcppExports.R"))

a <- commandArgs(trailingOnly = TRUE)
n  <- if (length(a) >= 1) as.integer(a[1]) else 256L
m  <- if (length(a) >= 2) as.integer(a[2]) else 50L
out <- if (length(a) >= 3) a[3] else "/tmp/work/scaling.csv"
which <- if (length(a) >= 4) a[4] else "all"   # all | direct | iter
do_direct <- (n <= 600L) && which %in% c("all","direct")   # cap the direct factorization (memory/time) at ~360k cells
do_blockcg <- (n <= 200L) && which %in% c("all","iter")  # block-CG (Jacobi) is a small-graph method; skip above ~40k cells

nr <- n; nc <- n; N <- as.double(nr) * nc
# vectorized 8-neighbour edge list (1-based)
cells <- seq_len(N); row <- ((cells - 1L) %/% nc) + 1L; col <- ((cells - 1L) %% nc) + 1L
from <- c(cells[col < nc], cells[row < nr], cells[row < nr & col < nc], cells[row < nr & col > 1])
to   <- c(cells[col < nc] + 1L, cells[row < nr] + nc, cells[row < nr & col < nc] + nc + 1L,
          cells[row < nr & col > 1] + nc - 1L)
ep <- cbind(as.integer(from), as.integer(to)); storage.mode(ep) <- "integer"
set.seed(1)
# heterogeneous conductance: smooth gradient + a low-conductance barrier column
barrier <- as.numeric(col == (nc %/% 2L))
cond <- as.numeric(exp(0.4 * scale(row) - 2.0 * barrier))
demes <- as.integer(round(seq(1, N - 1, length.out = m)))
Z <- matrix(-1 / N, N - 1L, m); kd <- demes < N
Z[cbind(demes[kd], which(kd))] <- Z[cbind(demes[kd], which(kd))] + 1

tt <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- force(expr); list(t = proc.time()[["elapsed"]] - t0, v = v) }

res <- list(n = n, N = N, edges = nrow(ep), m = m)
Xd <- NULL; fill <- NA
if (do_direct) {
  L <- forceSymmetric(assemble_reduced_laplacian(cond, ep))
  fac <- tt(Matrix::Cholesky(L))
  sol <- tt(as.matrix(Matrix::solve(fac$v, Z)))
  Xd <- sol$v
  fill <- length(as(fac$v, "Matrix")@x) / length(L@x)
  res$direct_setup_s <- round(fac$t, 3); res$direct_solve_s <- round(sol$t, 3)
  res$fill_ratio <- round(fill, 2)
} else { res$direct_setup_s <- NA; res$direct_solve_s <- NA; res$fill_ratio <- NA }

if (which %in% c("all","iter")) {
  ic <- tt(pcg_reduced_laplacian_ic(Z, cond, ep, tol = 1e-6, maxit = 5000))
  res$icpcg_s <- round(ic$t, 3); res$icpcg_iters <- max(ic$v$iterations)
  res$icpcg_conv <- all(ic$v$converged)
  res$icpcg_err  <- if (!is.null(Xd)) signif(max(abs(ic$v$solution - Xd)), 3) else NA
} else { res$icpcg_s <- NA; res$icpcg_iters <- NA; res$icpcg_conv <- NA; res$icpcg_err <- NA }

if (do_blockcg) {
  bc <- tt(block_cg_reduced_laplacian(Z, cond, ep, tol = 1e-6, maxit = 5000))
  res$blockcg_s <- round(bc$t, 3); res$blockcg_iters <- bc$v$iterations
  res$blockcg_conv <- all(bc$v$converged)
  res$blockcg_err  <- if (!is.null(Xd)) signif(max(abs(bc$v$solution - Xd)), 3) else NA
} else { res$blockcg_s <- NA; res$blockcg_iters <- NA; res$blockcg_conv <- NA; res$blockcg_err <- NA }

row <- as.data.frame(res, stringsAsFactors = FALSE)
write.table(row, out, sep = ",", append = file.exists(out),
            col.names = !file.exists(out), row.names = FALSE)
print(row)
