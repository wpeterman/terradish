# Scenario 5 validation: with deme sizes set so Ne ~ pi^k (effective size genuinely
# tied to the migration-stationary distribution), does the FRAME-coupled model
# recover BOTH the directional coefficient (~BETA) and the coupling exponent
# (alpha ~ k)? Contrast with uniform (ignores drift) and decoupled drift.
# Run from the terradish repo root:  Rscript dev/slim/validate_scen5.R
suppressMessages(devtools::load_all(".", quiet = TRUE)); suppressMessages(library(terra))

OUT <- Sys.getenv("OUTDIR", "dev/slim/out"); PRE <- file.path(OUT, "scen5_ts")
Y <- as.matrix(read.csv(paste0(PRE, "_Y.csv"), header = FALSE))
N <- scan(paste0(PRE, "_N.csv"), quiet = TRUE)
demeids <- scan(paste0(PRE, "_demeids.csv"), quiet = TRUE)
demes <- read.csv(paste0(PRE, "_demes.csv"))
tj <- readLines(paste0(PRE, "_truth.json"))
beta_true <- as.numeric(sub('.*"beta": *([-0-9.]+).*', "\\1", grep('"beta"', tj, value = TRUE)[1]))
k_true    <- as.numeric(sub('.*"k": *([-0-9.]+).*',    "\\1", grep('"k"',    tj, value = TRUE)[1]))

MAF <- 0.05
fr <- colSums(Y) / sum(N); Y <- Y[, pmin(fr, 1 - fr) >= MAF, drop = FALSE]
S <- cov_from_biallelic(Y, N = N); nu <- ncol(Y)

rd <- match(demeids, demes$deme)
coords <- cbind(demes$gx[rd], demes$gy[rd]); el <- demes$elev[rd]
DIM <- length(unique(demes$gx))
r <- rast(nrows = DIM, ncols = DIM, xmin = min(demes$gx) - .5, xmax = max(demes$gx) + .5,
          ymin = min(demes$gy) - .5, ymax = max(demes$gy) + .5)
vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- el; elevr <- setValues(r, vals)
graph <- conductance_surface(elevr, coords, directions = 4); vc <- graph$vertex_coordinates
perm <- apply(vc, 1, function(xy) which(abs(coords[, 1] - xy[1]) < 1e-6 & abs(coords[, 2] - xy[2]) < 1e-6))
dirn <- vc[, 1]; S_g <- S[perm, perm]
cat(sprintf("scen5: %d demes, %d markers; truth gamma_dir=%.2f  alpha(=k)=%.2f\n",
            nrow(vc), nu, beta_true, k_true))

fu <- dragon(dirn, graph, S_g, nu, "uniform")
fd <- dragon(dirn, graph, S_g, nu, "drift", drift = dirn)     # decoupled drift, covariate=elev
fc <- dragon(dirn, graph, S_g, nu, "coupled")
cat(sprintf("  uniform : gamma_dir=%+.3f\n", coef(fu)["gamma_dir"]))
cat(sprintf("  drift   : gamma_dir=%+.3f  kappa=%+.3f\n", coef(fd)["gamma_dir"], coef(fd)["kappa"]))
cat(sprintf("  coupled : gamma_dir=%+.3f  alpha=%+.3f   <- expect gamma~%.2f, alpha~%.2f\n",
            coef(fc)["gamma_dir"], coef(fc)["alpha"], beta_true, k_true))
cat(sprintf("  AIC  uniform=%.1f  drift=%.1f  coupled=%.1f\n", AIC(fu), AIC(fd), AIC(fc)))
cat("\nReading: scen5 sets Ne ~ pi^k, so the coupled model is correctly specified.\n",
    "If coupled recovers gamma_dir~BETA AND alpha in (0,1) near k (not pinned at the\n",
    "bound as in scen4), the end-to-end fair test passes: coupling works when Ne\n",
    "genuinely tracks the migration-stationary distribution.\n")
