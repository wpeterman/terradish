# Validate the packaged directed engine: numDeriv on the engine gradient,
# reduction to symmetric at gamma=0, and terradish_directed() recovery.
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra)); suppressMessages(library(numDeriv))
set.seed(1)

data(melip)
alt <- terra::unwrap(melip.altitude); fc <- terra::unwrap(melip.forestcover)
co  <- terra::unwrap(melip.coords)
covs <- c(terra::scale(alt), terra::scale(fc)); names(covs) <- c("altitude", "forestcover")
# coarsen so the per-absorber R solves are fast (engine math already proven;
# this just confirms the packaging). Full-resolution fitting needs the C++ backend.
covs <- terra::aggregate(covs, fact = 10, fun = "mean", na.rm = TRUE)
keep <- 1:10
surface <- conductance_surface(covs, co[keep, ], directions = 8)
nf <- length(surface$demes)

# directional covariate from (coarse, scaled) altitude
dir_cov <- edge_gradient(covs[["altitude"]], surface)
gen <- terradish:::.directed_generator(~ forestcover, surface, dir_cov)
cat(sprintf("nodes=%d, focal=%d, directed edges=%d, p=%d q=%d\n",
            gen$n, nf, nrow(gen$edges), gen$p, gen$q))

# simulate a Wishart distance response from the directed model
theta_true <- 0.5; gamma_true <- 0.6; nu <- 800
E <- terradish_directed_algorithm(gen, NULL, surface, S = NULL,
                                  par = c(theta_true, gamma_true))$covariance
Sigma <- 1.0 * E + 0.2 * diag(nf)
set.seed(7); Scov <- rWishart(1, df = nu, Sigma = Sigma / nu)[, , 1]
S <- outer(diag(Scov), rep(1, nf)) + outer(rep(1, nf), diag(Scov)) - 2 * Scov  # distance
diag(S) <- 0

cat("\n=== (0) ISOLATION: forward+adjoint only, L = sum(W*E), no measurement model ===\n")
# This tests the directed forward map and the transpose-solve adjoint in
# isolation. dL/dE = W (fixed symmetric); analytic grad uses the same internal
# adjoint via a tiny measurement-model-free path.
set.seed(3); W <- matrix(rnorm(nf^2), nf, nf); W <- (W + t(W)) / 2
adj_only <- function(par) {
  fwd <- terradish_directed_algorithm(gen, NULL, surface, S = NULL, par = par)
  # reproduce the adjoint chain with dL/dE = W
  Jc <- diag(nf) - matrix(1/nf, nf, nf)
  dL_dR <- -0.5 * (Jc %*% W %*% Jc); dL_dH <- dL_dR + t(dL_dR)
  ed <- gen$edges; ea <- ed[,1]; eb <- ed[,2]; n <- gen$n; focal <- surface$demes
  dL_drate <- numeric(length(ea))
  for (fj in seq_len(nf)) {
    hc <- fwd$hcache[[fj]]; j <- hc$j; idx <- hc$idx
    b <- numeric(n); b[focal] <- dL_dH[, fj]; b[j] <- 0
    ared <- as.numeric(Matrix::solve(Matrix::t(hc$Q), b[idx]))
    afull <- numeric(n); afull[idx] <- ared
    keep <- ea != j
    dL_drate[keep] <- dL_drate[keep] + afull[ea[keep]] * (hc$hfull[ea[keep]] - hc$hfull[eb[keep]])
  }
  rate <- gen$rates(par)
  c(as.numeric(crossprod(gen$sab, dL_drate * rate)), as.numeric(crossprod(gen$D, dL_drate * rate)))
}
L_of_par <- function(par) sum(W * terradish_directed_algorithm(gen, NULL, surface, S = NULL, par = par)$covariance)
par <- c(0.3, 0.4)
an0 <- adj_only(par); nd0 <- numDeriv::grad(L_of_par, par)
cat(sprintf("analytic: %s\nnumeric : %s\nmax abs diff: %.2e  %s\n",
            paste(sprintf("%.5f", an0), collapse = " "),
            paste(sprintf("%.5f", nd0), collapse = " "),
            max(abs(an0 - nd0)), ifelse(max(abs(an0 - nd0)) < 1e-4, "OK", "**FAIL**")))

cat("\n=== numDeriv check via generalized_wishart (full engine gradient) ===\n")
an <- terradish_directed_algorithm(gen, generalized_wishart, surface, S, par, nu = nu, gradient = TRUE)
nd <- numDeriv::grad(function(pp)
  terradish_directed_algorithm(gen, generalized_wishart, surface, S, pp, nu = nu, gradient = FALSE)$objective, par)
cat(sprintf("analytic: %s\nnumeric : %s\nmax abs diff: %.2e  %s\n",
            paste(sprintf("%.5f", an$gradient), collapse = " "),
            paste(sprintf("%.5f", nd), collapse = " "),
            max(abs(an$gradient - nd)), ifelse(max(abs(an$gradient - nd)) < 1e-4, "OK", "**FAIL**")))

cat("\n=== gamma=0 => reversible generator (symmetric rates) ===\n")
m <- nrow(gen$edges) / 2
r0 <- gen$rates(c(0.4, 0))
cat(sprintf("max|rate(a->b) - rate(b->a)| at gamma=0 = %.2e (should be ~0)\n",
            max(abs(r0[seq_len(m)] - r0[m + seq_len(m)]))))

cat("\n=== terradish_directed() recovery ===\n")
fit <- terradish_directed(S ~ forestcover, data = surface, directional = dir_cov,
                          measurement_model = generalized_wishart, nu = nu)
print(fit)
cat(sprintf("\ntheta_hat=%.3f (true %.2f) ; gamma_hat=%.3f (true %.2f)\n",
            fit$theta[1], theta_true, fit$gamma[1], gamma_true))
cat("\nDONE\n")
