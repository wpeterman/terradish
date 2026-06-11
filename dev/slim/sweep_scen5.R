# scen5 bias/variance sweep over the true directional coefficient (BETA) and the
# coupling exponent (k, = alpha). For each (BETA, k, seed) cell the recap driver
# simulates Ne ~ pi^k end-to-end; we fit uniform + coupled and record recovered
# gamma_dir and alpha. terradish is loaded ONCE; only the recap is shelled per cell.
# Writes dev/slim/scen5_sweep.csv incrementally. Run from the terradish repo root.
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(terra))

PY     <- Sys.getenv("PY", "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe")
RECAP  <- "dev/slim/run_scen5_recap.py"
OUTDIR <- "dev/slim/sweep_tmp"
RESCSV <- Sys.getenv("RESCSV", "dev/slim/scen5_sweep.csv")
NGEN   <- Sys.getenv("NGEN", "10000"); K0 <- Sys.getenv("K0", "300")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

build_graph <- function(demes, demeids) {
  rd <- match(demeids, demes$deme)
  coords <- cbind(demes$gx[rd], demes$gy[rd]); el <- demes$elev[rd]
  DIM <- length(unique(demes$gx))
  r <- rast(nrows = DIM, ncols = DIM, xmin = min(demes$gx) - .5, xmax = max(demes$gx) + .5,
            ymin = min(demes$gy) - .5, ymax = max(demes$gy) + .5)
  vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- el
  graph <- conductance_surface(setValues(r, vals), coords, directions = 4)
  vc <- graph$vertex_coordinates
  perm <- apply(vc, 1, function(xy) which(abs(coords[, 1] - xy[1]) < 1e-6 & abs(coords[, 2] - xy[2]) < 1e-6))
  list(graph = graph, perm = perm, dirn = vc[, 1])
}

fit_cell <- function(beta, k, seed) {
  out <- system2(PY, c(RECAP, "--beta", beta, "--k", k, "--seed", seed,
                       "--out", OUTDIR, "--ngen", NGEN, "--k0", K0),
                 stdout = TRUE, stderr = TRUE)
  PRE <- file.path(OUTDIR, "scen5_ts")
  if (!file.exists(paste0(PRE, "_Y.csv"))) stop("recap produced no output")
  Y <- as.matrix(read.csv(paste0(PRE, "_Y.csv"), header = FALSE))
  N <- scan(paste0(PRE, "_N.csv"), quiet = TRUE)
  demeids <- scan(paste0(PRE, "_demeids.csv"), quiet = TRUE)
  demes <- read.csv(paste0(PRE, "_demes.csv"))
  fr <- colSums(Y) / sum(N); Y <- Y[, pmin(fr, 1 - fr) >= 0.05, drop = FALSE]
  S <- cov_from_biallelic(Y, N = N); nu <- ncol(Y)
  g <- build_graph(demes, demeids); S_g <- S[g$perm, g$perm]
  fu <- dragon(g$dirn, g$graph, S_g, nu, "uniform", n_start = 8L)
  fc <- dragon(g$dirn, g$graph, S_g, nu, "coupled", n_start = 8L)
  data.frame(beta = beta, k = k, seed = seed, nu = nu,
             gdir_uniform = unname(coef(fu)["gamma_dir"]),
             gdir_coupled = unname(coef(fc)["gamma_dir"]),
             alpha_coupled = unname(coef(fc)["alpha"]))
}

grid <- expand.grid(beta = 0.5, k = c(0.0, 0.25, 0.5, 0.75), seed = 1:4)
cat(sprintf("scen5 sweep: %d cells (NGEN=%s K0=%s)\n", nrow(grid), NGEN, K0))
for (i in seq_len(nrow(grid))) {
  cell <- grid[i, ]
  row <- tryCatch(fit_cell(cell$beta, cell$k, cell$seed),
                  error = function(e) { cat("  cell", i, "failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(row)) next
  write.table(row, RESCSV, sep = ",", row.names = FALSE, col.names = (i == 1),
              append = (i != 1))
  cat(sprintf("  [%d/%d] beta=%.1f k=%.2f seed=%d -> gdir_u=%+.3f gdir_c=%+.3f alpha=%+.3f\n",
              i, nrow(grid), cell$beta, cell$k, cell$seed,
              row$gdir_uniform, row$gdir_coupled, row$alpha_coupled))
}

cat("\n===== summary (mean +/- sd over seeds) =====\n")
res <- read.csv(RESCSV)
agg <- aggregate(cbind(gdir_uniform, gdir_coupled, alpha_coupled) ~ beta + k, res,
                 function(x) c(m = mean(x), s = sd(x)))
for (i in seq_len(nrow(agg))) {
  cat(sprintf("beta=%.1f k=%.2f : gdir_c=%+.2f+/-%.2f  alpha=%+.2f+/-%.2f  (alpha interior? %s)\n",
              agg$beta[i], agg$k[i], agg$gdir_coupled[i, "m"], agg$gdir_coupled[i, "s"],
              agg$alpha_coupled[i, "m"], agg$alpha_coupled[i, "s"],
              ifelse(agg$alpha_coupled[i, "m"] > 0.02 && agg$alpha_coupled[i, "m"] < 0.98, "yes", "NO")))
}
cat("\nwrote", RESCSV, "\n")
