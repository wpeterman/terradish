# Diagnose Scenario 3 null result: does the symmetric engine find IBR (tau>0) on
# the same data, and is the directed engine's gradient nonzero away from (0,0)?
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))

OUTPRE <- "dev/slim/scen3ts"; DIM <- 5L
Y <- as.matrix(read.csv(paste0(OUTPRE, "_Y.csv"), header = FALSE))
N <- scan(paste0(OUTPRE, "_N.csv"), quiet = TRUE)
demeids <- scan(paste0(OUTPRE, "_demeids.csv"), quiet = TRUE)
demes <- read.csv(paste0(OUTPRE, "_demes.csv"))
elev_by_deme <- demes$elev[match(demeids, demes$deme)]
coords <- cbind(demes$gx, demes$gy)[match(demeids, demes$deme), ]
fr <- colSums(Y) / sum(N); Y <- Y[, pmin(fr, 1 - fr) >= 0.05, drop = FALSE]
S <- cov_from_biallelic(Y, N = N)
nu <- 100

r <- rast(nrows = DIM, ncols = DIM, xmin = -0.5, xmax = DIM - 0.5, ymin = -0.5, ymax = DIM - 0.5)
vals <- rep(NA_real_, ncell(r)); vals[cellFromXY(r, coords)] <- elev_by_deme
elevr <- setValues(r, vals); names(elevr) <- "elev"
surface <- conductance_surface(elevr, coords, directions = 8)
dir_cov <- edge_gradient(elevr, surface)

cat("=== (1) standard symmetric terradish on the same S ===\n")
fit_sym <- suppressWarnings(terradish(S ~ elev, data = surface,
                      conductance_model = loglinear_conductance,
                      measurement_model = wishart_covariance, nu = nu, leverage = FALSE,
                      control = NewtonRaphsonControl(maxit = 100, verbose = FALSE)))
cat(sprintf("symmetric: theta_elev=%.3f  phi(tau,sigma)=[%s]  boundary=%s  loglik=%.2f\n",
            coef(fit_sym)[1], paste(round(fit_sym$fit$phi[, 1], 3), collapse = ", "),
            fit_sym$fit$boundary, fit_sym$loglik))

cat("\n=== (2) directed engine: phi (tau) and gradient away from (0,0) ===\n")
gen <- terradish:::.directed_generator(~ elev, surface, dir_cov)
for (par in list(c(0, 0), c(0.3, 0.3), c(0.5, 0.5))) {
  r2 <- terradish_directed_algorithm(gen, wishart_covariance, surface, S, par = par,
                                     nu = nu, gradient = TRUE)
  cat(sprintf("par=(%.1f,%.1f): obj=%.2f  phi=[%s]  boundary=%s  grad=[%s]\n",
              par[1], par[2], r2$objective,
              paste(round(r2$phi, 3), collapse = ","), r2$boundary,
              paste(round(r2$gradient, 3), collapse = ",")))
}

cat("\n=== (3) does commute-time E (gamma=0) correlate with S? ===\n")
E0 <- terradish_directed_algorithm(gen, NULL, surface, NULL, par = c(0, 0))$covariance
lt <- lower.tri(S)
cat(sprintf("cor(E_commute_offdiag, S_offdiag) = %.3f\n", cor(E0[lt], S[lt])))
cat(sprintf("cor(diag(E_commute), diag(S)) = %.3f\n", cor(diag(E0), diag(S))))
cat("DONE\n")
