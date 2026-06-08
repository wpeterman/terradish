# =====================================================================
# Tier 1 SLiM validation, Scenario 1 (stepping-stone, varying deme size).
# Runs SLiM, builds the genetic covariance, and fits the terradish drift
# surface vs a scalar nugget. Validates that the diagonal drift term
# recovers the Ne gradient (NOT the conductance/IBR term, since migration
# is uniform) and improves model fit.
#
# Self-contained: SLiM writes per-deme allele frequencies; finite sampling
# is applied in R (no Python / tree-sequence dependency).
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

SLIM   <- "C:/msys64/mingw64/bin/slim.exe"
SCRIPT <- "dev/slim/scen1_stepping_stone.slim"
OUTPRE <- "dev/slim/scen1"
DIM    <- 5
NSAMP  <- 25
SEED   <- 1001
run_slim <- !file.exists(paste0(OUTPRE, "_freq.csv")) ||
            isTRUE(as.logical(Sys.getenv("SLIM_FORCE", "FALSE")))

if (run_slim) {
  cat("Running SLiM scenario 1 ...\n")
  # Leaner than a full SNP panel so the forward-only run completes quickly
  # (no recapitation available here); convertToSubstitution keeps the mutation
  # registry small. Drift gradient is strong (BETA=0.7 => ~4x deme-size contrast).
  args <- c("-seed", SEED,
            "-d", paste0("DIM=", DIM),
            "-d", "N0=180", "-d", "BETA=0.7", "-d", "MIG=0.03",
            "-d", "MU=1e-7", "-d", "RHO=1e-8", "-d", "L=400000",
            "-d", paste0("NSAMP=", NSAMP), "-d", "BURNIN=2500",
            "-d", paste0("OUTPRE='", OUTPRE, "'"),
            SCRIPT)
  st <- system2(SLIM, args, stdout = TRUE, stderr = TRUE)
  cat(tail(st, 3), sep = "\n"); cat("\n")
}

## ---- read SLiM output ----
freq  <- read.csv(paste0(OUTPRE, "_freq.csv"))
demes <- read.csv(paste0(OUTPRE, "_demes.csv"))
ndeme <- nrow(demes)
P <- t(as.matrix(freq[, grep("^deme", colnames(freq))]))   # demes x loci (true freqs)
cat(sprintf("demes=%d, raw loci=%d\n", ndeme, ncol(P)))

## ---- finite sampling: Y ~ Binom(2*NSAMP, p) ----
set.seed(99)
nchrom <- 2 * NSAMP
Y <- matrix(rbinom(length(P), nchrom, P), nrow = nrow(P))   # demes x loci
# keep loci polymorphic in the pooled sample
fr <- colSums(Y) / (nchrom * ndeme)
keep <- fr > 0 & fr < 1
Y <- Y[, keep, drop = FALSE]
nloc <- ncol(Y)
cat(sprintf("polymorphic sampled loci=%d\n", nloc))

S <- cov_from_biallelic(Y, N = nchrom)

## ---- ground-truth diagnostic: does within-site variance track Ne? ----
cat(sprintf("diag(S) range: [%.2f, %.2f]; mean |off-diag|: %.3f\n",
            min(diag(S)), max(diag(S)), mean(abs(S[lower.tri(S)]))))
cat(sprintf("cor(diag(S), cov)  = %.3f  (expected NEGATIVE: more cov -> more Ne -> less drift)\n",
            cor(diag(S), demes$cov)))
cat(sprintf("cor(diag(S), size) = %.3f  (expected NEGATIVE)\n",
            cor(diag(S), demes$size)))

## ---- terradish surface on the DIM x DIM grid ----
cov_by_deme <- demes$cov                       # drift covariate, deme order
coords <- cbind(demes$gx, demes$gy)
r <- rast(nrows = DIM, ncols = DIM,
          xmin = -0.5, xmax = DIM - 0.5, ymin = -0.5, ymax = DIM - 0.5)
vals <- rep(NA_real_, ncell(r))
cells <- cellFromXY(r, coords)
vals[cells] <- cov_by_deme
gradient <- setValues(r, vals)
names(gradient) <- "cov"
surface <- conductance_surface(gradient, coords, directions = 8)

## ---- choose nu (independent-marker count); point estimates are nu-invariant ----
nu <- max(20, round(nloc / 20))   # rough LD-thinned independent count
cat(sprintf("nu (for SE/AIC only) = %d\n", nu))

## ---- fit drift surface vs scalar nugget ----
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
cat("\n================ RESULTS (Scenario 1) ================\n")
cat("True: deme size = 250*exp(0.7*cov)  ->  Ne increases with cov\n")
cat("      => drift/nugget should DECREASE with cov  => gamma_cov < 0\n\n")
cat("drift model phi:\n"); print(round(fit_d$fit$phi, 3))
cat(sprintf("\nconductance theta (IBR effect of cov; expect ~0): %.3f\n",
            fit_d$mle$theta["cov"]))
cat(sprintf("drift slope gamma_cov: %.3f  (expected NEGATIVE)\n",
            fit_d$fit$phi["gamma_cov", 1]))
cat(sprintf("AIC drift=%.1f  scalar=%.1f  -> drift preferred: %s\n",
            aic(fit_d), aic(fit_s), aic(fit_d) < aic(fit_s)))
cat("======================================================\n")
