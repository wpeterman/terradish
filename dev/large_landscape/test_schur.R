# Validation/unit test for the tiled Schur prototype (schur_tiled.R).
# Demonstrates EXACT focal-site resistance recovery (contrast: prototype B's
# node-lumping biased E by ~25%), for both single-shot and tiled elimination.
# Run: Rscript test_schur.R   (Matrix only)
suppressMessages(library(Matrix))
HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "schur_tiled.R"))

## build a 30x30 grid with a low-conductance vertical barrier (heterogeneous) ---
nr <- 30L; nc <- 30L; N <- nr * nc
cid <- function(r, c) (r - 1L) * nc + c
E <- list()
for (r in 1:nr) for (cc in 1:nc) {
  if (cc < nc) E[[length(E) + 1L]] <- c(cid(r, cc), cid(r, cc + 1L))
  if (r < nr)  E[[length(E) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc))
  if (r < nr && cc < nc) E[[length(E) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc + 1L))
  if (r < nr && cc > 1)  E[[length(E) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc - 1L))
}
edges <- do.call(rbind, E)
rr <- ((seq_len(N) - 1L) %/% nc) + 1L; ccs <- ((seq_len(N) - 1L) %% nc) + 1L
logc <- 0.4 * scale(rr) - 2.5 * (ccs == 15L)          # barrier at column 15
cond <- as.numeric(exp(logc - mean(logc)))
w <- cond[edges[, 1]] + cond[edges[, 2]]
L <- full_laplacian(N, edges, w)
focal <- as.integer(c(cid(3, 3), cid(5, 25), cid(20, 8), cid(27, 27), cid(15, 2), cid(12, 28)))

## reference: effective resistance among focal from the FULL graph -------------
R_full <- eff_resistance(L, focal, focal[1])

## single-shot Schur (eliminate all non-focal) --------------------------------
red <- schur_eliminate(L, setdiff(seq_len(N), focal))
foc_local <- match(focal, red$keep)
R_schur <- eff_resistance(Matrix(red$S, sparse = TRUE), foc_local, foc_local[1])
err_schur <- max(abs(R_schur - R_full)) / max(abs(R_full))

## tiled Schur (eliminate interior in 4 quadrant tiles) -----------------------
tile_of <- (rr > nr / 2) * 2L + (ccs > nc / 2) + 1L
tiles <- split(seq_len(N), tile_of)
tl <- schur_tiled(L, focal, tiles)
foc_local2 <- match(focal, tl$keep)
R_tiled <- eff_resistance(Matrix(tl$S, sparse = TRUE), foc_local2, foc_local2[1])
err_tiled <- max(abs(R_tiled - R_full)) / max(abs(R_full))
max_tile_interior <- max(sapply(tiles, function(t) length(setdiff(t, focal))))

cat(sprintf("full graph cells                 : %d\n", N))
cat(sprintf("focal sites                      : %d\n", length(focal)))
cat(sprintf("single-shot Schur rel error      : %.2e  (exact up to round-off)\n", err_schur))
cat(sprintf("tiled Schur rel error            : %.2e  (== single-shot; Schur composes)\n", err_tiled))
cat(sprintf("largest per-tile interior solve  : %d cells (vs %d full)\n", max_tile_interior, N - 1L))
cat(sprintf("contrast: prototype B node-lumping biased E by ~25%%; Schur is exact.\n"))

stopifnot(err_schur < 1e-8, err_tiled < 1e-8,
          max(abs(R_tiled - R_schur)) < 1e-8)
cat("PASS: focal resistance preserved exactly by single-shot and tiled elimination.\n")
