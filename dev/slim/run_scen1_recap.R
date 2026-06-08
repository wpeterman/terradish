# =====================================================================
# Scenario 1 via the RECAPITATION pipeline (tree-seq + pyslim + msprime).
# SLiM records a tree sequence; Python recapitates, overlays neutral mutations,
# and samples per deme; R builds the covariance and fits drift vs scalar.
# Deep coalescent ancestry removes the forward-only young-mutation artifact, so
# we expect a well-scaled covariance, a clean drift gradient, AND tau > 0.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

SLIM   <- "C:/msys64/mingw64/bin/slim.exe"
PY     <- "C:/Users/peterman.73/AppData/Local/anaconda3/envs/slim/python.exe"
SLIM_SCRIPT <- "dev/slim/scen1_ts.slim"
PY_SCRIPT   <- "dev/slim/recap_sample.py"
OUTPRE <- "dev/slim/scen1ts"
DIM    <- 5
NSAMP  <- 25
SEED   <- 1001
force  <- isTRUE(as.logical(Sys.getenv("SLIM_FORCE", "FALSE")))

if (force || !file.exists(paste0(OUTPRE, "_Y.csv"))) {
  cat("1) SLiM tree-seq ...\n")
  # Low migration -> detectable Fst (drift differences + genuine IBR => tau>0);
  # long burn-in (cheap with tree-seq) so within-deme coalescence completes.
  system2(SLIM, c("-seed", SEED, "-d", paste0("DIM=", DIM),
                  "-d", "N0=150", "-d", "BETA=0.9", "-d", "MIG=0.005",
                  "-d", "RHO=1e-8", "-d", "L=1000000", "-d", "BURNIN=4000",
                  "-d", paste0("OUTPRE='", OUTPRE, "'"), SLIM_SCRIPT),
          stdout = TRUE, stderr = TRUE) |> tail(2) |> cat(sep = "\n"); cat("\n")

  demes0 <- read.csv(paste0(OUTPRE, "_demes.csv"))
  anc_Ne <- sum(demes0$size)
  cat(sprintf("2) recapitate (ancestral_Ne=%d) + mutations + sample ...\n", anc_Ne))
  out <- system2(PY, c(PY_SCRIPT, "--trees", paste0(OUTPRE, ".trees"),
                       "--out", OUTPRE, "--ancestral_Ne", anc_Ne,
                       "--recombination_rate", "1e-8", "--mutation_rate", "1e-8",
                       "--nsample", NSAMP, "--seed", SEED, "--mode", "population"),
                 stdout = TRUE, stderr = TRUE)
  cat(tail(out, 2), sep = "\n"); cat("\n")
}

## ---- read recapitated genotype summaries ----
Y <- as.matrix(read.csv(paste0(OUTPRE, "_Y.csv"), header = FALSE))   # demes x loci
N <- scan(paste0(OUTPRE, "_N.csv"), quiet = TRUE)                    # per-deme haploid
demeids <- scan(paste0(OUTPRE, "_demeids.csv"), quiet = TRUE)
demes <- read.csv(paste0(OUTPRE, "_demes.csv"))
cov_by_deme <- demes$cov[match(demeids, demes$deme)]
coords <- cbind(demes$gx, demes$gy)[match(demeids, demes$deme), ]
ndeme <- nrow(Y)
cat(sprintf("demes=%d, loci=%d, N range=[%d,%d]\n", ndeme, ncol(Y), min(N), max(N)))

# SNP QC: minor-allele-frequency filter (recapitated data has a full SFS; rare
# variants otherwise dominate the standardized covariance). MAF >= 0.05 is a
# standard landscape-genetics threshold.
MAF <- as.numeric(Sys.getenv("MAF", "0.05"))
fr  <- colSums(Y) / sum(N)
maf <- pmin(fr, 1 - fr)
keep <- maf >= MAF
Y <- Y[, keep, drop = FALSE]
cat(sprintf("loci after MAF>=%.2f filter: %d (of %d)\n", MAF, ncol(Y), length(maf)))
S <- cov_from_biallelic(Y, N = N)

## ---- ground-truth diagnostics ----
cat(sprintf("diag(S) range: [%.2f, %.2f]; mean |off-diag|: %.3f\n",
            min(diag(S)), max(diag(S)), mean(abs(S[lower.tri(S)]))))
cat(sprintf("cor(diag(S), cov)  = %.3f  (expected NEGATIVE)\n", cor(diag(S), cov_by_deme)))

## ---- terradish surface + fit ----
r <- rast(nrows = DIM, ncols = DIM, xmin = -0.5, xmax = DIM - 0.5,
          ymin = -0.5, ymax = DIM - 0.5)
vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- cov_by_deme
gradient <- setValues(r, vals); names(gradient) <- "cov"
surface <- conductance_surface(gradient, coords, directions = 8)
nu <- max(20, round(ncol(Y) / 20))

g_drift <- wishart_drift_covariates(data.frame(cov = cov_by_deme),
                                    model = "wishart_covariance")
fit_d <- suppressWarnings(terradish(S ~ cov, data = surface,
                   conductance_model = loglinear_conductance,
                   measurement_model = g_drift, nu = nu, leverage = FALSE,
                   control = NewtonRaphsonControl(maxit = 200, verbose = FALSE)))
fit_s <- suppressWarnings(terradish(S ~ cov, data = surface,
                   conductance_model = loglinear_conductance,
                   measurement_model = wishart_covariance, nu = nu, leverage = FALSE,
                   control = NewtonRaphsonControl(maxit = 200, verbose = FALSE)))
aic <- function(f) -2 * f$loglik + 2 * (length(f$mle$theta) + nrow(f$fit$phi))

cat("\n========= RESULTS (Scenario 1, recapitated) =========\n")
cat("drift model phi:\n"); print(round(fit_d$fit$phi, 3))
cat(sprintf("tau (IBR scale): %.3f  (expect > 0 now)\n", fit_d$fit$phi["tau", 1]))
cat(sprintf("gamma_cov (drift slope): %.3f  (expect NEGATIVE)\n", fit_d$fit$phi["gamma_cov", 1]))
cat(sprintf("AIC drift=%.1f scalar=%.1f -> drift preferred: %s\n",
            aic(fit_d), aic(fit_s), aic(fit_d) < aic(fit_s)))
cat("=====================================================\n")
