# Localize the coarse-melip adjoint mismatch: run the PACKAGE directed engine on
# a small, well-conditioned raster (smooth covariates, rook adjacency). If the
# isolation test matches here, the package code is correct and the coarse-melip
# failure is numerical conditioning (wide covariate range -> exp rate blow-up).
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra)); suppressMessages(library(numDeriv))
set.seed(1)

mk <- function(DIM, scale_cov, directions) {
  r <- rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM, ymin = 0, ymax = DIM)
  gx <- xFromCell(r, seq_len(ncell(r))); gy <- yFromCell(r, seq_len(ncell(r)))
  v1 <- scale_cov * scale(gx + 0.5 * gy)[, 1]   # smooth, controlled range
  elev <- scale_cov * scale(gx)[, 1]
  covs <- c(setValues(r, v1), setValues(r, elev)); names(covs) <- c("v1", "elev")
  fc <- c(1L, DIM, DIM * DIM, DIM * (DIM - 1) + 1L, (DIM * DIM) %/% 2L)
  coords <- xyFromCell(r, fc)
  surface <- conductance_surface(covs, coords, directions = directions)
  dir_cov <- edge_gradient(covs[["elev"]], surface)
  gen <- terradish:::.directed_generator(~ v1, surface, dir_cov)
  list(surface = surface, gen = gen, nf = length(surface$demes))
}

isolation <- function(o, par) {
  gen <- o$gen; surface <- o$surface; nf <- o$nf
  set.seed(3); W <- matrix(rnorm(nf^2), nf, nf); W <- (W + t(W)) / 2
  Lf <- function(p) sum(W * terradish_directed_algorithm(gen, NULL, surface, NULL, p)$covariance)
  fwd <- terradish_directed_algorithm(gen, NULL, surface, NULL, par)
  Jc <- diag(nf) - matrix(1/nf, nf, nf)
  dL_dR <- -0.5 * (Jc %*% W %*% Jc); dL_dH <- dL_dR + t(dL_dR)
  ed <- gen$edges; ea <- ed[,1]; eb <- ed[,2]; n <- gen$n; focal <- surface$demes
  dr <- numeric(length(ea))
  for (fj in seq_len(nf)) {
    hc <- fwd$hcache[[fj]]; j <- hc$j; idx <- hc$idx
    b <- numeric(n); b[focal] <- dL_dH[, fj]; b[j] <- 0
    ared <- as.numeric(Matrix::solve(Matrix::t(hc$Q), b[idx]))
    afull <- numeric(n); afull[idx] <- ared
    keep <- ea != j
    dr[keep] <- dr[keep] + afull[ea[keep]] * (hc$hfull[ea[keep]] - hc$hfull[eb[keep]])
  }
  rate <- gen$rates(par)
  an <- c(as.numeric(crossprod(gen$sab, dr * rate)), as.numeric(crossprod(gen$D, dr * rate)))
  nd <- numDeriv::grad(Lf, par)
  list(an = an, nd = nd, diff = max(abs(an - nd)))
}

for (cfg in list(list(DIM=6, sc=1.0, dir=4, par=c(0.3,0.4)),
                 list(DIM=6, sc=1.0, dir=8, par=c(0.3,0.4)),
                 list(DIM=8, sc=1.0, dir=8, par=c(0.3,0.4)),
                 list(DIM=8, sc=2.5, dir=8, par=c(0.3,0.4)))) {
  o <- mk(cfg$DIM, cfg$sc, cfg$dir)
  res <- isolation(o, cfg$par)
  cat(sprintf("DIM=%d scale=%.1f dir=%d nodes=%d nf=%d: max|an-nd|=%.2e  %s\n",
              cfg$DIM, cfg$sc, cfg$dir, o$gen$n, o$nf, res$diff,
              ifelse(res$diff < 1e-4, "OK", "**FAIL**")))
}
cat("DONE\n")
