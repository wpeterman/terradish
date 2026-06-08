# =====================================================================
# Scenario 2 via the RECAPITATION pipeline (continuous space).
# SLiM (tree-seq) -> pyslim.recapitate -> msprime mutations -> spatial Voronoi
# sampling in Python -> covariance -> fit drift vs scalar. Tests whether the
# drift surface recovers a habitat-driven local-Ne gradient in continuous space
# with realistic (deep-coalescent) ancestry.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

SLIM   <- "C:/msys64/mingw64/bin/slim.exe"
PY     <- "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe"
SLIM_SCRIPT <- "dev/slim/scen2_ts.slim"
PY_SCRIPT   <- "dev/slim/recap_sample.py"
OUTPRE <- "dev/slim/scen2ts"
HAB0     <- 0.15
HABSLOPE <- 0.85
K        <- 2800
NSAMP    <- 30
SEED     <- 2002
force  <- isTRUE(as.logical(Sys.getenv("SLIM_FORCE", "FALSE")))

if (force || !file.exists(paste0(OUTPRE, "_Y.csv"))) {
  cat("1) SLiM tree-seq (continuous space) ...\n")
  system2(SLIM, c("-seed", SEED, "-d", paste0("K=", K),
                  "-d", "SIGMA_C=0.10", "-d", "SIGMA_M=0.10", "-d", "SIGMA_D=0.08",
                  "-d", "RHO=1e-8", "-d", "L=1000000", "-d", "FGRID=5",
                  "-d", "BURNIN=4000",
                  "-d", paste0("HAB0=", HAB0), "-d", paste0("HABSLOPE=", HABSLOPE),
                  "-d", paste0("OUTPRE='", OUTPRE, "'"), SLIM_SCRIPT),
          stdout = TRUE, stderr = TRUE) |> tail(2) |> cat(sep = "\n"); cat("\n")

  anc_Ne <- round(K * (HAB0 + HABSLOPE * 0.5))
  cat(sprintf("2) recapitate (ancestral_Ne=%d) + mutations + spatial sample ...\n", anc_Ne))
  out <- system2(PY, c(PY_SCRIPT, "--trees", paste0(OUTPRE, ".trees"),
                       "--out", OUTPRE, "--ancestral_Ne", anc_Ne,
                       "--recombination_rate", "1e-8", "--mutation_rate", "1e-8",
                       "--nsample", NSAMP, "--seed", SEED, "--mode", "spatial",
                       "--focal", paste0(OUTPRE, "_focal.csv")),
                 stdout = TRUE, stderr = TRUE)
  cat(tail(out, 2), sep = "\n"); cat("\n")
}

Y <- as.matrix(read.csv(paste0(OUTPRE, "_Y.csv"), header = FALSE))
N <- scan(paste0(OUTPRE, "_N.csv"), quiet = TRUE)
demeids <- scan(paste0(OUTPRE, "_demeids.csv"), quiet = TRUE)
demes <- read.csv(paste0(OUTPRE, "_demes.csv"))
hab_by_deme <- demes$habitat[match(demeids, demes$deme)]
coords <- cbind(demes$fx, demes$fy)[match(demeids, demes$deme), ]
cat(sprintf("focal demes=%d, loci=%d, N range=[%d,%d]\n", nrow(Y), ncol(Y), min(N), max(N)))

MAF <- as.numeric(Sys.getenv("MAF", "0.05"))
fr  <- colSums(Y) / sum(N)
maf <- pmin(fr, 1 - fr)
Y <- Y[, maf >= MAF, drop = FALSE]
cat(sprintf("loci after MAF>=%.2f: %d\n", MAF, ncol(Y)))
S <- cov_from_biallelic(Y, N = N)

cat(sprintf("diag(S) range: [%.2f, %.2f]; mean |off-diag|: %.3f\n",
            min(diag(S)), max(diag(S)), mean(abs(S[lower.tri(S)]))))
cat(sprintf("cor(diag(S), habitat) = %.3f  (expected NEGATIVE)\n",
            cor(diag(S), hab_by_deme)))

r <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
xc <- xFromCell(r, seq_len(ncell(r)))
habitat <- setValues(r, HAB0 + HABSLOPE * xc); names(habitat) <- "habitat"
surface <- conductance_surface(habitat, coords, directions = 8)
nu <- max(20, round(ncol(Y) / 20))

g_drift <- wishart_drift_covariates(data.frame(habitat = hab_by_deme),
                                    model = "wishart_covariance")
fit_d <- suppressWarnings(terradish(S ~ habitat, data = surface,
                   conductance_model = loglinear_conductance,
                   measurement_model = g_drift, nu = nu, leverage = FALSE,
                   control = NewtonRaphsonControl(maxit = 200, verbose = FALSE)))
fit_s <- suppressWarnings(terradish(S ~ habitat, data = surface,
                   conductance_model = loglinear_conductance,
                   measurement_model = wishart_covariance, nu = nu, leverage = FALSE,
                   control = NewtonRaphsonControl(maxit = 200, verbose = FALSE)))
aic <- function(f) -2 * f$loglik + 2 * (length(f$mle$theta) + nrow(f$fit$phi))

cat("\n========= RESULTS (Scenario 2, continuous space, recapitated) =========\n")
cat("drift model phi:\n"); print(round(fit_d$fit$phi, 3))
cat(sprintf("tau (IBR scale): %.3f\n", fit_d$fit$phi["tau", 1]))
cat(sprintf("gamma_habitat (drift slope): %.3f  (expect NEGATIVE)\n",
            fit_d$fit$phi["gamma_habitat", 1]))
cat(sprintf("AIC drift=%.1f scalar=%.1f -> drift preferred: %s\n",
            aic(fit_d), aic(fit_s), aic(fit_d) < aic(fit_s)))
cat("=======================================================================\n")
