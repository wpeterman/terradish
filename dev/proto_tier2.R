# =====================================================================
# Tier 2 prototype: hierarchical conductance field  log c = X theta + Z_u u
# Validates the core mechanism before hardening into package code:
#  (1) penalized gradient over (theta, u) is correct (numDeriv),
#  (2) tau^2 -> small reduces to plain terradish,
#  (3) recovery: an UNMAPPED conductance feature is absorbed by the field u
#      while the mapped covariate effects theta are protected.
# Optimizer: L-BFGS over (theta, u) (gradient-only) -> avoids the dense Hessian.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(terra))
suppressMessages(library(numDeriv))
set.seed(20)

## ---- conductance model directly from a design matrix (loglinear) ----
loglinear_from_matrix <- function(D) {
  D <- as.matrix(D)
  cm <- function(theta) {
    cond <- as.vector(exp(D %*% theta))
    dfm  <- cond * D
    list(conductance = cond,
         confint = function(...) NULL,
         df__dx = function(k) cond * theta[k],
         df__dtheta = function(k) dfm[, k],
         df__dtheta_matrix = dfm,
         d2f__dtheta_dtheta = function(k, l) cond * D[, k] * D[, l],
         d2f__dtheta_dx = function(k, l) cond * ((k == l) + D[, k] * theta[l]))
  }
  class(cm) <- c("terradish_conductance_model", "radish_conductance_model")
  attr(cm, "default") <- rep(0, ncol(D))
  cm
}

## ---- coarse piecewise-constant field basis + GMRF precision ----
build_field <- function(coords, G = 6, eps = 1e-3) {
  xr <- range(coords[, 1]); yr <- range(coords[, 2])
  # tiny pad so max falls inside the last bin
  cx <- pmin(G, 1 + floor(G * (coords[, 1] - xr[1]) / (diff(xr) + 1e-9)))
  cy <- pmin(G, 1 + floor(G * (coords[, 2] - yr[1]) / (diff(yr) + 1e-9)))
  cell <- (cy - 1) * G + cx                       # 1..G*G
  occ  <- sort(unique(cell))
  m    <- length(occ)
  remap <- match(cell, occ)
  Z <- matrix(0, nrow(coords), m)
  Z[cbind(seq_len(nrow(coords)), remap)] <- 1
  # coarse-grid rook adjacency among occupied cells
  gx <- ((occ - 1) %% G) + 1; gy <- ((occ - 1) %/% G) + 1
  A <- matrix(0, m, m)
  for (a in 1:m) for (b in 1:m)
    if (a != b && abs(gx[a] - gx[b]) + abs(gy[a] - gy[b]) == 1) A[a, b] <- 1
  Lc <- diag(rowSums(A)) - A
  list(Z = Z, Q = Lc + eps * diag(m), m = m, occ = occ, G = G,
       cell_xy = cbind(gx, gy))
}

## ---- melip surface ----
data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)
covs <- c(terra::scale(melip.altitude), terra::scale(melip.forestcover))
names(covs) <- c("altitude", "forestcover")
surface <- conductance_surface(covs, melip.coords, directions = 8)
n  <- length(surface$demes)
Vc <- surface$vertex_coordinates
X  <- as.matrix(surface$x[, c("altitude", "forestcover")])
N  <- nrow(X)

## ---- TRUE surface: mapped covariates + an UNMAPPED smooth blob ----
xr <- range(Vc[, 1]); yr <- range(Vc[, 2])
bx <- xr[1] + 0.35 * diff(xr); by <- yr[1] + 0.6 * diff(yr)
rad <- 0.18 * sqrt(diff(xr)^2 + diff(yr)^2)
blob <- 1.6 * exp(-(((Vc[, 1] - bx)^2 + (Vc[, 2] - by)^2) / (2 * rad^2)))
blob <- blob - mean(blob)                        # unmapped log-conductance feature
theta_true <- c(altitude = 0.5, forestcover = -0.6)
D_true <- cbind(X, blob = blob)
cm_true <- loglinear_from_matrix(D_true)
E_true <- as.matrix(terradish_algorithm(
  cm_true, leastsquares, surface, S = diag(n),
  theta = c(theta_true, 1.0), objective = FALSE,
  gradient = FALSE, hessian = FALSE, partial = FALSE)$covariance)

tau_true <- 1.0; sigma_true <- 0.25; nu <- 1500
S <- rWishart(1, df = nu, Sigma = (tau_true * E_true + sigma_true * diag(n)) / nu)[, , 1]

## ---- field + penalized objective over (theta, u) ----
fld <- build_field(Vc, G = 6)
cat(sprintf("field cells m = %d\n", fld$m))
D_fit <- cbind(X, fld$Z)
cm_fit <- loglinear_from_matrix(D_fit)
p <- ncol(X); uidx <- (p + 1):(p + fld$m)

pen_fn <- function(par, tau2, grad = FALSE) {
  res <- terradish_algorithm(cm_fit, wishart_covariance, surface, S = S,
                             theta = par, nu = nu, gradient = grad,
                             hessian = FALSE, partial = FALSE)
  u <- par[uidx]
  obj <- res$objective + as.numeric(t(u) %*% fld$Q %*% u) / (2 * tau2)
  if (!grad) return(obj)
  g <- res$gradient
  g[uidx] <- g[uidx] + as.vector(fld$Q %*% u) / tau2
  list(obj = obj, grad = g)
}

## ---- (1) numDeriv check of the penalized gradient ----
set.seed(5)
par0 <- c(0.2, -0.3, rnorm(fld$m, 0, 0.2))
tau2 <- 0.5
an <- pen_fn(par0, tau2, grad = TRUE)
nd <- numDeriv::grad(function(p) pen_fn(p, tau2, grad = FALSE), par0)
cat(sprintf("(1) penalized gradient max|analytic-numeric| = %.2e  %s\n",
            max(abs(an$grad - nd)),
            ifelse(max(abs(an$grad - nd)) < 1e-4, "OK", "**FAIL**")))

## ---- (3) fit by L-BFGS at a fixed tau^2 ----
fit_lbfgs <- function(tau2) {
  optim(par0, fn = function(p) pen_fn(p, tau2, grad = FALSE),
        gr = function(p) pen_fn(p, tau2, grad = TRUE)$grad,
        method = "L-BFGS-B",
        control = list(maxit = 300, factr = 1e7))
}
opt <- fit_lbfgs(0.5)
theta_hat <- opt$par[1:p]; u_hat <- opt$par[uidx]
u_cells <- as.vector(fld$Z %*% u_hat)            # field mapped back to cells
cat(sprintf("\n(3) recovery at tau2=0.5 (convergence code %d):\n", opt$convergence))
cat(sprintf("    theta_alt  = %.3f (true %.2f)\n", theta_hat[1], theta_true[1]))
cat(sprintf("    theta_fc   = %.3f (true %.2f)\n", theta_hat[2], theta_true[2]))
cat(sprintf("    cor(field, unmapped blob) = %.3f  (want strongly positive)\n",
            cor(u_cells, blob)))

## ---- compare: no-field fit (omitted-variable bias baseline) ----
cm_nofield <- loglinear_from_matrix(X)
fit_nf <- terradish(S ~ altitude + forestcover, data = surface,
                    conductance_model = loglinear_conductance,
                    measurement_model = wishart_covariance, nu = nu,
                    leverage = FALSE,
                    control = NewtonRaphsonControl(maxit = 100, verbose = FALSE))
cat(sprintf("\n    no-field theta_alt = %.3f, theta_fc = %.3f (true %.2f, %.2f)\n",
            coef(fit_nf)[1], coef(fit_nf)[2], theta_true[1], theta_true[2]))

## ---- (2) tau^2 -> 0 shrinks the field ----
opt_small <- fit_lbfgs(1e-4)
cat(sprintf("\n(2) tau2=1e-4: max|u| = %.4f (should be ~0, reducing to plain terradish)\n",
            max(abs(opt_small$par[uidx]))))
cat("\nDONE\n")
