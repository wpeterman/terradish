# =====================================================================
# Scenario 3 (Tier 3 validation): asymmetric migration via the recapitation
# pipeline. SLiM (tree-seq, downstream-biased migration) -> recapitate ->
# mutations -> per-deme covariance -> fit the DIRECTED model (full: theta+gamma)
# vs the REDUCED reversible model (gamma=0) on the same commute-time engine.
# Expect gamma recovered (sign +) and the directed model preferred (LRT/AIC) -
# the empirical test of whether direction is identifiable from real coalescent
# data. Runs on a small lattice graph where the R directed engine is fast.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

SLIM <- "C:/msys64/mingw64/bin/slim.exe"
PY   <- "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe"
OUTPRE <- "dev/slim/scen3ts"; DIM <- 5L; N0 <- 200L; NSAMP <- 25L; SEED <- 3003L
BETA_DIR <- 0.8
force <- isTRUE(as.logical(Sys.getenv("SLIM_FORCE", "FALSE")))

if (force || !file.exists(paste0(OUTPRE, "_Y.csv"))) {
  cat("1) SLiM tree-seq (asymmetric migration) ...\n")
  system2(SLIM, c("-seed", SEED, "-d", paste0("DIM=", DIM), "-d", paste0("N0=", N0),
                  "-d", "MIG=0.004", "-d", paste0("BETA_DIR=", BETA_DIR),
                  "-d", "RHO=1e-8", "-d", "L=1000000", "-d", "BURNIN=6000",
                  "-d", paste0("OUTPRE='", OUTPRE, "'"), "dev/slim/scen3_ts.slim"),
          stdout = TRUE, stderr = TRUE) |> tail(2) |> cat(sep = "\n"); cat("\n")
  anc <- DIM * DIM * N0
  cat(sprintf("2) recapitate (ancestral_Ne=%d) + mutations + sample ...\n", anc))
  system2(PY, c("dev/slim/recap_sample.py", "--trees", paste0(OUTPRE, ".trees"),
                "--out", OUTPRE, "--ancestral_Ne", anc, "--recombination_rate", "1e-8",
                "--mutation_rate", "1e-8", "--nsample", NSAMP, "--seed", SEED,
                "--mode", "population"),
          stdout = TRUE, stderr = TRUE) |> tail(2) |> cat(sep = "\n"); cat("\n")
}

Y <- as.matrix(read.csv(paste0(OUTPRE, "_Y.csv"), header = FALSE))
N <- scan(paste0(OUTPRE, "_N.csv"), quiet = TRUE)
demeids <- scan(paste0(OUTPRE, "_demeids.csv"), quiet = TRUE)
demes <- read.csv(paste0(OUTPRE, "_demes.csv"))
elev_by_deme <- demes$elev[match(demeids, demes$deme)]
coords <- cbind(demes$gx, demes$gy)[match(demeids, demes$deme), ]

MAF <- 0.05
fr <- colSums(Y) / sum(N); Y <- Y[, pmin(fr, 1 - fr) >= MAF, drop = FALSE]
S <- cov_from_biallelic(Y, N = N)
cat(sprintf("demes=%d, loci(MAF>=%.2f)=%d\n", nrow(Y), MAF, ncol(Y)))
cat(sprintf("diag(S) range [%.2f, %.2f]; mean|off-diag| %.3f\n",
            min(diag(S)), max(diag(S)), mean(abs(S[lower.tri(S)]))))

## terradish surface (small lattice) + elevation directional covariate
r <- rast(nrows = DIM, ncols = DIM, xmin = -0.5, xmax = DIM - 0.5, ymin = -0.5, ymax = DIM - 0.5)
vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- elev_by_deme
elevr <- setValues(r, vals); names(elevr) <- "elev"
surface <- conductance_surface(elevr, coords, directions = 8)
dir_cov <- edge_gradient(elevr, surface)
nu <- max(20, round(ncol(Y) / 20))

## full directed fit (theta on elevation [symmetric] + gamma [directional])
fit_dir <- terradish_directed(S ~ elev, data = surface, directional = dir_cov,
                              measurement_model = wishart_covariance, nu = nu)
## reduced reversible fit: gamma fixed at 0 (same engine)
gen <- terradish:::.directed_generator(~ elev, surface, dir_cov)
ll_reduced <- function() {
  f <- function(th) terradish_directed_algorithm(gen, wishart_covariance, surface, S,
                                                 par = c(th, 0), nu = nu, gradient = FALSE)$objective
  o <- optimize(f, c(-5, 5)); -o$objective
}
llr <- ll_reduced()
llf <- fit_dir$loglik
lrt <- 2 * (llf - llr); pval <- pchisq(lrt, df = 1, lower.tail = FALSE)

cat("\n========= RESULTS (Scenario 3: asymmetric migration) =========\n")
print(fit_dir)
cat(sprintf("\ntau (IBR scale) = %.3f  (0 => no spatial signal detected)\n",
            fit_dir$phi["tau"]))
cat(sprintf("gamma_hat (directional) = %.3f  (expect POSITIVE: downhill-biased)\n", fit_dir$gamma[1]))
cat(sprintf("LRT directed vs reversible: 2dLL=%.2f  df=1  p=%.4g  (small p => direction detected)\n",
            lrt, pval))
cat("==============================================================\n")
