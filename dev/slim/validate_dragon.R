# =====================================================================
# Validate dragon() against a REAL conductance_surface() graph, reusing the
# scen4 demes. This exercises the one piece the reference tests cannot: the
# wrapper's alignment of S, the directional potential, and the drift covariate
# with the graph's active-cell (vertex) ordering.
#
# Run from the terradish repo root:
#   Rscript dev/slim/validate_dragon.R            # default: scen4 corr, seed 1
#   VARIANT=uncorr OUTDIR=dev/slim/out Rscript dev/slim/validate_dragon.R
# Requires the scen4 outputs (run dev/slim/run_scen4_recap.py first), terra, and
# the DRAGON engine in R/dragon.R.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

VARIANT <- Sys.getenv("VARIANT", "corr")                 # "corr" (collinear) or "uncorr"
OUTDIR  <- Sys.getenv("OUTDIR",  "dev/slim/out")          # seed folder from run_scen4_recap.py
PRE     <- file.path(OUTDIR, paste0("scen4_", VARIANT, "_ts"))
MAF     <- 0.05
stopifnot(file.exists(paste0(PRE, "_Y.csv")))

## ---- 1. genetic covariance from the recap output (S is in demeids order) ----
Y <- as.matrix(read.csv(paste0(PRE, "_Y.csv"), header = FALSE))
N <- scan(paste0(PRE, "_N.csv"), quiet = TRUE)
demeids <- scan(paste0(PRE, "_demeids.csv"), quiet = TRUE)        # sampled deme ids
demes   <- read.csv(paste0(PRE, "_demes.csv"))                    # deme, gx, gy, elev, dens

fr <- colSums(Y) / sum(N)
Y  <- Y[, pmin(fr, 1 - fr) >= MAF, drop = FALSE]
S  <- cov_from_biallelic(Y, N = N)                               # rows/cols in demeids order
nu <- ncol(Y)

## per-S-row deme attributes (elevation = directional potential; dens = drift)
row_deme <- match(demeids, demes$deme)
S_gx   <- demes$gx[row_deme];  S_gy <- demes$gy[row_deme]
S_elev <- demes$elev[row_deme]; S_dens <- demes$dens[row_deme]

## ---- 2. build a real conductance_surface graph on the 5x5 deme lattice ------
DIM <- length(unique(demes$gx))
r <- rast(nrows = DIM, ncols = DIM,
          xmin = min(demes$gx) - 0.5, xmax = max(demes$gx) + 0.5,
          ymin = min(demes$gy) - 0.5, ymax = max(demes$gy) + 0.5)
coords <- cbind(S_gx, S_gy)
vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- S_elev
elevr <- setValues(r, vals); names(elevr) <- "elev"
graph <- conductance_surface(elevr, coords, directions = 4)      # rook adjacency

## ---- 3. ALIGN S / covariates to the graph's vertex order (the key step) -----
vc <- graph$vertex_coordinates                                   # graph-vertex coords
# perm[v] = the S row whose deme sits at graph vertex v
perm <- apply(vc, 1, function(xy)
  which(abs(coords[, 1] - xy[1]) < 1e-6 & abs(coords[, 2] - xy[2]) < 1e-6))
stopifnot(length(perm) == nrow(vc), !anyNA(perm), !any(duplicated(perm)))

S_g     <- S[perm, perm]
directional <- vc[, 1]                                           # x-coord = elevation
drift_g <- S_dens[perm]
# sanity: the directional potential must equal the matched demes' elevation
stopifnot(max(abs(directional - S_elev[perm])) < 1e-6)
cat(sprintf("graph: %d demes, %d edges; markers (MAF>=%.2f) = %d\n",
            nrow(vc), nrow(graph$edge_pairs), MAF, nu))

## ---- 4. fit DRAGON three ways + model comparison ----------------------------
fit_blind <- dragon(directional, graph, S_g, nu, coalescence = "uniform")
fit_drift <- dragon(directional, graph, S_g, nu, coalescence = "drift", drift = drift_g)
fit_coupl <- dragon(directional, graph, S_g, nu, coalescence = "coupled")

lrt_drift <- 2 * (fit_drift$loglik - fit_blind$loglik)          # drift vs blind, df=1

cat("\n================ DRAGON on a real graph: ", VARIANT, " ================\n", sep = "")
cat("drift-blind :"); print(round(coef(fit_blind), 3))
cat("drift-joint :"); print(round(coef(fit_drift), 3))
cat("coupled     :"); print(round(coef(fit_coupl), 3))
cat(sprintf("\nLRT drift vs blind = %.2f (df=1; >3.84 => drift detected)\n", lrt_drift))
cat(sprintf("AIC  blind=%.1f  drift=%.1f  coupled=%.1f\n",
            AIC(fit_blind), AIC(fit_drift), AIC(fit_coupl)))

## ---- 5. collinearity diagnostic --------------------------------------------
diag <- dragon_collinearity(directional, graph, S_g, nu, drift = drift_g)
cat(sprintf("\ncollinearity: r_cov=%.2f  gamma_dir blind=%.3f joint=%.3f gap=%.2f\n",
            diag$r_cov, diag$gdir_blind, diag$gdir_joint, diag$gap))
cat("  ->", diag$verdict, "\n")

## ---- 6. what to expect (cross-check vs the Python scen4 analysis) -----------
cat("\nExpected (see DRAGON/dev/slim/SCEN4_REVIEW.md):\n")
if (VARIANT == "corr") {
  cat("  collinear drift: r_cov ~ 1.0, diagnostic should WARN CONFOUNDED;\n",
      "  drift-blind gamma_dir positive (direction recovered, ~0.7-1.6 across seeds);\n",
      "  drift-joint pulls gamma_dir DOWN (drift steals directional signal).\n")
} else {
  cat("  orthogonal drift: r_cov ~ 0, diagnostic OK; drift-joint ~ drift-blind\n",
      "  (gamma_dir stable near the blind estimate).\n")
}
cat("\nIf the directional sign, the collinearity verdict, and the blind-vs-joint\n",
    "pattern match the Python analysis, the dragon() wrapper is aligned correctly.\n")

## ---- 7. (optional) contrast with the commute-time engine on the same graph --
## Uncomment to show DRAGON (coalescent) vs terradish_directed (commute time):
# dir_cov <- edge_gradient(elevr, graph)
# fit_commute <- terradish_directed(S_g ~ elev, data = graph, directional = dir_cov,
#                                   measurement_model = wishart_covariance, nu = nu)
# cat("\ncommute-time gamma:", round(fit_commute$gamma[1], 3),
#     " (expect attenuated vs DRAGON's gamma_dir)\n")
