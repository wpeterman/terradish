# =====================================================================
# Tier 3, Phase 1-2 prototype: directed-generator engine + reverse-mode adjoint.
# Validates the forward map (theta,gamma -> commute-time covariance E) and the
# analytic gradient dL/d(theta,gamma) against numDeriv, before hardening into
# package code. Crux of Tier 3: the transpose-solve adjoint.
# =====================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(Matrix))
suppressMessages(library(numDeriv))
set.seed(1)

## ---- small graph: 5x5 lattice, focal subset, node covariate + directional cov ----
DIM <- 5L; n <- DIM * DIM
gx <- rep(0:(DIM - 1L), times = DIM); gy <- rep(0:(DIM - 1L), each = DIM)
coords <- cbind(gx, gy)
X <- as.numeric(scale(gx + 0.5 * gy))          # node covariate (symmetric part)
elev <- as.numeric(gx)                          # elevation -> directional covariate
focal <- c(1L, 5L, 11L, 13L, 21L, 25L)          # focal node indices
nf <- length(focal)

ee <- which(as.matrix(dist(coords)) == 1, arr.ind = TRUE)
ee <- ee[ee[, 1] != ee[, 2], , drop = FALSE]    # directed edges (both ways)
ea <- ee[, 1]; eb <- ee[, 2]
dab <- elev[ea] - elev[eb]                       # directional covariate (antisymmetric)
sab <- (X[ea] + X[eb]) / 2                       # symmetric endpoint covariate

## ---- forward: (theta, gamma) -> E (commute-time covariance among focal) ----
Jc <- diag(nf) - matrix(1 / nf, nf, nf)

build_G <- function(theta, gamma) {
  rate <- exp(theta * sab + gamma * dab)
  G <- sparseMatrix(i = ea, j = eb, x = rate, dims = c(n, n))
  diag(G) <- 0
  G <- G - Diagonal(n, x = rowSums(G))
  list(G = G, rate = rate)
}

forward <- function(theta, gamma) {
  gen <- build_G(theta, gamma)
  G <- gen$G
  H <- matrix(0, n, nf)            # hitting times: node -> focal absorber
  hcache <- vector("list", nf)
  for (fj in seq_len(nf)) {
    j <- focal[fj]
    idx <- setdiff(seq_len(n), j)
    Q <- G[idx, idx, drop = FALSE]
    hred <- tryCatch(as.numeric(solve(Q, rep(-1, n - 1L))),
                     error = function(e) NULL)
    if (is.null(hred)) return(NULL)             # ill-conditioned (strong asymmetry)
    hfull <- numeric(n); hfull[idx] <- hred
    H[, fj] <- hfull
    hcache[[fj]] <- list(j = j, idx = idx, Q = Q, hfull = hfull)
  }
  Hf <- H[focal, , drop = FALSE]   # nf x nf : H[i,j]=h(focal_i -> focal_j)
  R <- Hf + t(Hf)                  # symmetric commute time
  E <- -0.5 * Jc %*% R %*% Jc
  list(E = as.matrix(E), gen = gen, hcache = hcache, Hf = Hf, R = R)
}

## ---- reverse-mode adjoint: dL/dE -> dL/d(theta,gamma) ----
adjoint <- function(fwd, dL_dE) {
  dL_dR <- -0.5 * (Jc %*% dL_dE %*% Jc)          # E = -1/2 J R J
  dL_dH <- dL_dR + t(dL_dR)                       # R = Hf + Hf^T
  dL_drate <- numeric(length(ea))
  for (fj in seq_len(nf)) {
    hc <- fwd$hcache[[fj]]; j <- hc$j; idx <- hc$idx
    # b over transient nodes: sensitivities of L to H[focal_i, fj]
    b <- numeric(n)
    b[focal] <- dL_dH[, fj]
    b[j] <- 0
    bred <- b[idx]
    ared <- as.numeric(solve(t(hc$Q), bred))      # TRANSPOSE solve
    afull <- numeric(n); afull[idx] <- ared
    hfull <- hc$hfull
    # dL/drate(a->b) += a[a]*(h[a]-h[b]) for edges with a != j
    keep <- ea != j
    dL_drate[keep] <- dL_drate[keep] +
      afull[ea[keep]] * (hfull[ea[keep]] - hfull[eb[keep]])
  }
  # chain to (theta, gamma): d rate/d theta = rate*sab ; d rate/d gamma = rate*dab
  rate <- fwd$gen$rate
  dL_dtheta <- sum(dL_drate * rate * sab)
  dL_dgamma <- sum(dL_drate * rate * dab)
  c(theta = dL_dtheta, gamma = dL_dgamma)
}

## ---- simulate data from the model, then check forward + adjoint ----
theta_true <- 0.4; gamma_true <- 0.5
Etrue <- forward(theta_true, gamma_true)$E
nu <- 300
# distance response for leastsquares: S = double-centering inverse of E (a distance)
Rtrue <- outer(diag(Etrue), rep(1, nf)) + outer(rep(1, nf), diag(Etrue)) - 2 * Etrue
set.seed(7); S <- Rtrue + matrix(rnorm(nf * nf, 0, 0.05 * sd(Rtrue)), nf, nf)
S <- (S + t(S)) / 2; diag(S) <- 0

profile <- function(par, want_grad = FALSE) {
  fwd <- forward(par[1], par[2])
  if (is.null(fwd)) {                            # ill-conditioned region
    if (!want_grad) return(1e10)
    return(list(obj = 1e10, grad = c(theta = 0, gamma = 0)))
  }
  sub <- radish_subproblem(leastsquares, fwd$E, S, nu = nu,
                           control = NewtonRaphsonControl(verbose = FALSE,
                                                          ftol = 1e-10, ctol = 1e-10))
  if (!want_grad) return(sub$loglikelihood)
  list(obj = sub$loglikelihood, grad = adjoint(fwd, as.matrix(sub$gradient)))
}

cat("=== forward sanity ===\n")
cat(sprintf("E symmetric: %s ; PSD-ish min eig = %.3g\n",
            isTRUE(all.equal(Etrue, t(Etrue))),
            min(eigen(Etrue, symmetric = TRUE, only.values = TRUE)$values)))

cat("\n=== adjoint vs numDeriv ===\n")
par0 <- c(0.25, 0.3)
an <- profile(par0, want_grad = TRUE)
nd <- numDeriv::grad(function(p) profile(p, FALSE), par0)
cat(sprintf("analytic grad: theta=%.5f gamma=%.5f\n", an$grad[1], an$grad[2]))
cat(sprintf("numeric  grad: theta=%.5f gamma=%.5f\n", nd[1], nd[2]))
cat(sprintf("max abs diff: %.3e  %s\n", max(abs(an$grad - nd)),
            ifelse(max(abs(an$grad - nd)) < 1e-4, "OK", "**FAIL**")))

cat("\n=== gamma=0 reduces to a reversible (symmetric) generator ===\n")
g0 <- build_G(0.3, 0)$G
asym <- max(abs(g0 - t(g0)))
cat(sprintf("at gamma=0, max|G - G^T| = %.3e (0 => reversible/symmetric rates)\n", asym))

cat("\n=== recovery: fit (theta,gamma) by bounded optim ===\n")
opt <- optim(c(0, 0), function(p) profile(p, FALSE),
             gr = function(p) profile(p, TRUE)$grad, method = "L-BFGS-B",
             lower = c(-2, -2), upper = c(2, 2))
cat(sprintf("theta_hat=%.3f (true %.2f) ; gamma_hat=%.3f (true %.2f) ; conv=%d\n",
            opt$par[1], theta_true, opt$par[2], gamma_true, opt$convergence))

cat("\nDONE\n")
