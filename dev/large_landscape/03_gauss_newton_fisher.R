# =============================================================================
# Prototype D: Gauss-Newton / Fisher-information curvature (replaces exact Hessian)
# -----------------------------------------------------------------------------
# Attacks W3 (the (1+p)-solve exact Hessian) AND is the inference object.
#
# Idea: for the Gaussian/Wishart/least-squares measurement models the objective
# is a (transform of a) Gaussian negative log-likelihood, so the EXPECTED Hessian
# is the Fisher information
#       I(theta) = J^T W J ,   J = d vec(E) / d theta  (q x p),  W = m.m. weight
# where q = number of modelled E-entries. J's columns come from the SAME adjoint
# machinery as the gradient (one solve each), or tr-type terms can be estimated
# by Hutchinson probing (see hutchinson_diag.cpp). I(theta) is:
#   * positive semidefinite -> always a descent direction (more robust than the
#     possibly-indefinite exact Hessian far from the optimum);
#   * the asymptotic information matrix -> vcov = I^{-1}, SEs = sqrt(diag(vcov));
#   * cheaper: no per-parameter Hessian solve loop.
#
# This prototype shows, on a Gaussian model E ~ loss, that:
#   (i) the Gauss-Newton curvature J^T W J approximates the true Hessian and
#       equals it at the optimum;
#  (ii) Hutchinson trace estimation recovers tr(B) for curvature probing.
#
# INTEGRATION: replaces the hessian branch in terradish_algorithm(); SEs feed
# vcov()/confint()/summary(). Keep the exact Hessian available for small-N checks.
#
# Depends: Matrix; sources 01_block_solver_adjoint.R.
# =============================================================================

suppressMessages(library(Matrix))
if (!exists("forward_E")) source("01_block_solver_adjoint.R")

## ---- Jacobian J = d vec(E) / d theta via the adjoint engine -----------------
# One linearised solve per parameter (could be replaced by probing). Here we use
# forward-mode over theta for clarity; production should reuse adjoint solves.
jacobian_E <- function(theta, edges, demes, cond_fun, dcond_fun, h = 1e-6) {
  p <- length(theta)
  E0 <- as.numeric(forward_E(theta, edges, demes, cond_fun)$E)
  J <- matrix(0, length(E0), p)
  for (k in seq_len(p)) {
    tp <- theta; tp[k] <- tp[k] + h
    tm <- theta; tm[k] <- tm[k] - h
    Ep <- as.numeric(forward_E(tp, edges, demes, cond_fun)$E)
    Em <- as.numeric(forward_E(tm, edges, demes, cond_fun)$E)
    J[, k] <- (Ep - Em) / (2 * h)
  }
  list(E0 = E0, J = J)
}

## ---- Gauss-Newton / Fisher information --------------------------------------
# W: weight matrix on vec(E) implied by the measurement model. For an iid
# Gaussian residual model W = (1/var) I; for Wishart it is the metric from the
# Sigma parameterisation (tau E + exp(sigma) I). nu scales the whole thing
# (effective sample size) -> SEs scale as 1/sqrt(nu), matching generalized_wishart.
fisher_information <- function(J, W, nu = 1) nu * crossprod(J, W %*% J)

standard_errors <- function(I) {
  vcov <- tryCatch(solve(I), error = function(e) MASS::ginv(I))
  sqrt(pmax(diag(vcov), 0))
}

## ---- Hutchinson trace primitive (curvature probing) -------------------------
# Estimates tr(B) for an implicit SPD operator given as a function apply_B(z).
# Use for tr(L^{-1} M) type terms without forming inverses. See verify.py (~1%
# at 2000 probes for tr(L^{-1})). For curvature you usually need few probes.
hutchinson_trace <- function(apply_B, n, n_probe = 200L) {
  est <- numeric(n_probe)
  for (t in seq_len(n_probe)) {
    z <- sample(c(-1, 1), n, replace = TRUE)
    est[t] <- sum(z * apply_B(z))
  }
  c(mean = mean(est), se = sd(est) / sqrt(n_probe))
}

## ============================ self-test ======================================
if (sys.nframe() == 0L) {
  set.seed(3)
  nr <- 8L; nc <- 8L; N <- nr * nc
  cid <- function(r, c) (r - 1L) * nc + c
  E_ <- list()
  for (r in 1:nr) for (cc in 1:nc) {
    if (cc < nc) E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r, cc + 1L))
    if (r < nr)  E_[[length(E_) + 1L]] <- c(cid(r, cc), cid(r + 1L, cc))
  }
  edges <- do.call(rbind, E_)
  demes <- c(1L, 20L, 45L, 64L); m <- length(demes)
  X <- matrix(rnorm(N * 2), N, 2)
  cond_fun  <- function(th) as.numeric(exp(X %*% th))
  dcond_fun <- function(th) { cc <- cond_fun(th); cc * X }
  theta <- c(0.2, -0.3)

  # Gaussian model: vec(E) ~ N(vec(target), s2 I); loss = 0.5/s2 * ||E-target||^2
  s2 <- 0.5
  set.seed(99); target <- as.numeric(crossprod(matrix(rnorm(m * m), m, m)))
  W <- diag(1 / s2, length(target))

  jc <- jacobian_E(theta, edges, demes, cond_fun, dcond_fun)
  I_gn <- fisher_information(jc$J, W, nu = 10)
  se   <- standard_errors(I_gn)

  # true Hessian of the loss by finite differences (for comparison)
  loss <- function(th) {
    e <- as.numeric(forward_E(th, edges, demes, cond_fun)$E)
    (0.5 / s2) * sum((e - target)^2)
  }
  H <- matrix(0, 2, 2); h <- 1e-4
  for (a in 1:2) for (b in 1:2) {
    tpp<-theta;tpp[a]<-tpp[a]+h;tpp[b]<-tpp[b]+h
    tpm<-theta;tpm[a]<-tpm[a]+h;tpm[b]<-tpm[b]-h
    tmp<-theta;tmp[a]<-tmp[a]-h;tmp[b]<-tmp[b]+h
    tmm<-theta;tmm[a]<-tmm[a]-h;tmm[b]<-tmm[b]-h
    H[a,b]<-(loss(tpp)-loss(tpm)-loss(tmp)+loss(tmm))/(4*h*h)
  }
  # Gauss-Newton curvature at nu=1 (no nu scaling) for like-for-like comparison
  I_gn1 <- fisher_information(jc$J, W, nu = 1)

  cat("Gauss-Newton curvature (nu=1):\n"); print(round(I_gn1, 3))
  cat("True Hessian (finite diff):\n");    print(round(H, 3))
  cat(sprintf("relative Frobenius diff: %.3f  (small => GN ~ Hessian; gap is the\n",
              norm(I_gn1 - H, "F") / norm(H, "F")))
  cat("  curvature from residual*second-derivative, ~0 near a good fit)\n")
  cat(sprintf("SEs at nu=10: %s\n", paste(sprintf('%.4f', se), collapse = " ")))

  # Hutchinson check: tr(L^{-1}) on this graph
  L <- build_reduced_laplacian(cond_fun(theta), edges)
  Lf <- Matrix::Cholesky(L)
  apply_Linv <- function(z) as.numeric(Matrix::solve(Lf, z))
  ht <- hutchinson_trace(apply_Linv, nrow(L), n_probe = 1500L)
  cat(sprintf("Hutchinson tr(L^-1): %.3f +/- %.3f  | true %.3f\n",
              ht["mean"], ht["se"], sum(diag(solve(as.matrix(L))))))
}
