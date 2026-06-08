# =====================================================================
# Tier 1 SLiM validation, Scenario 2 (continuous space, habitat-driven density).
# Runs SLiM, builds the genetic covariance among focal sampling points, and
# fits the terradish drift surface vs a scalar nugget. Validates that the
# diagonal drift term recovers the habitat (Ne) gradient in a realistic
# continuous-space setting.
#
# Self-contained: SLiM writes per-focal allele frequencies; finite sampling
# (with per-focal sample sizes) is applied in R. No Python dependency.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

SLIM   <- "C:/msys64/mingw64/bin/slim.exe"
SCRIPT <- "dev/slim/scen2_continuous_space.slim"
OUTPRE <- "dev/slim/scen2"
HAB0     <- 0.15
HABSLOPE <- 0.85
SEED   <- 2002
run_slim <- !file.exists(paste0(OUTPRE, "_freq.csv")) ||
            isTRUE(as.logical(Sys.getenv("SLIM_FORCE", "FALSE")))

if (run_slim) {
  cat("Running SLiM scenario 2 (continuous space) ...\n")
  args <- c("-seed", SEED,
            "-d", "K=2800", "-d", "SIGMA_C=0.10", "-d", "SIGMA_M=0.12",
            "-d", "SIGMA_D=0.12", "-d", "MU=1e-7", "-d", "RHO=1e-8",
            "-d", "L=1000000", "-d", "FGRID=5", "-d", "NSAMP=30",
            "-d", "BURNIN=3000",
            "-d", paste0("HAB0=", HAB0), "-d", paste0("HABSLOPE=", HABSLOPE),
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
cat(sprintf("focal points=%d, raw loci=%d, nsamp range=[%d,%d]\n",
            ndeme, ncol(P), min(demes$nsamp), max(demes$nsamp)))

## ---- finite sampling with per-focal sample sizes ----
set.seed(99)
nchrom <- 2 * demes$nsamp                          # per-deme haploid sample size
Y <- matrix(0L, nrow = ndeme, ncol = ncol(P))
for (d in seq_len(ndeme))
  Y[d, ] <- rbinom(ncol(P), nchrom[d], P[d, ])
fr <- colSums(Y) / sum(nchrom)
keep <- fr > 0 & fr < 1
Y <- Y[, keep, drop = FALSE]
nloc <- ncol(Y)
cat(sprintf("polymorphic sampled loci=%d\n", nloc))

S <- cov_from_biallelic(Y, N = nchrom)             # N is per-population (length ndeme)

## ---- ground-truth diagnostic: does within-site variance track habitat? ----
# drift (diagonal of S) should DECREASE with habitat (higher Ne -> less drift)
cat(sprintf("diag(S) range: [%.2f, %.2f]; mean |off-diag|: %.3f\n",
            min(diag(S)), max(diag(S)), mean(abs(S[lower.tri(S)]))))
cat(sprintf("cor(diag(S), habitat) = %.3f  (expected NEGATIVE)\n",
            cor(diag(S), demes$habitat)))
cat(sprintf("cor(diag(S), nsamp)   = %.3f  (sampling-size confound check)\n",
            cor(diag(S), demes$nsamp)))

## ---- terradish surface: raster with the habitat covariate over [0,1]^2 ----
hab_cov <- demes$habitat
coords  <- cbind(demes$fx, demes$fy)
r <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
xc <- xFromCell(r, seq_len(ncell(r)))
habitat <- setValues(r, HAB0 + HABSLOPE * xc)
names(habitat) <- "habitat"
surface <- conductance_surface(habitat, coords, directions = 8)

## ---- nu (independent-marker count); point estimates are nu-invariant ----
nu <- max(20, round(nloc / 20))
cat(sprintf("nu (for SE/AIC only) = %d\n", nu))

## ---- fit drift surface vs scalar nugget ----
g_drift <- wishart_drift_covariates(data.frame(habitat = hab_cov),
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
cat("\n================ RESULTS (Scenario 2, continuous space) ================\n")
cat("True: local density = K*habitat, habitat increases with x => Ne increases with x\n")
cat("      => drift/nugget should DECREASE with habitat => gamma_habitat < 0\n\n")
cat("drift model phi:\n"); print(round(fit_d$fit$phi, 3))
cat(sprintf("\nconductance theta (habitat effect on movement): %.3f\n",
            fit_d$mle$theta["habitat"]))
cat(sprintf("drift slope gamma_habitat: %.3f  (expected NEGATIVE)\n",
            fit_d$fit$phi["gamma_habitat", 1]))
cat(sprintf("AIC drift=%.1f  scalar=%.1f  -> drift preferred: %s\n",
            aic(fit_d), aic(fit_s), aic(fit_d) < aic(fit_s)))
cat("========================================================================\n")
