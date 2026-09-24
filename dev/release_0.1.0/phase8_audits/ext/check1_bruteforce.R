# Brute-force likelihood checks for the measurement-model extensions.
suppressMessages({library(terradish); library(Matrix)})
set.seed(42)
n <- 7
X <- matrix(rnorm(n * n), n, n)
E <- crossprod(X) / n + diag(n) * 0.3
z <- rnorm(n)
nu <- 30

# hand-coded Wishart log density for W ~ W_p(nu, V)
ldwish <- function(W, nu, V) {
  p <- nrow(W)
  0.5 * (nu - p - 1) * as.numeric(determinant(W)$modulus) -
    0.5 * sum(diag(solve(V, W))) - 0.5 * nu * p * log(2) -
    0.5 * nu * as.numeric(determinant(V)$modulus) -
    (p * (p - 1) / 4 * log(pi) + sum(lgamma((nu + 1 - seq_len(p)) / 2)))
}

cat("==== (a) wishart_covariates, covariance model ====\n")
g <- wishart_covariates(data.frame(env = z), model = "wishart_covariance")
K <- attr(g, "kernel_covariates")[, , 1]
zc <- z  # scale = FALSE is the shared pairwise-covariate default.
H <- diag(n) - 1/n
L <- cbind(diag(n - 1), -1)
cat("kernel == centered absdiff:", all.equal(K, -0.5 * H %*% abs(outer(zc, zc, "-")) %*% H, check.attributes = FALSE), "\n")
cat("kernel rank:", qr(K)$rank, "\n")
stopifnot(isTRUE(all.equal(K, -0.5 * H %*% abs(outer(zc, zc, "-")) %*% H, check.attributes = FALSE)))
Sigma0 <- 0.8 * E + 0.5 * K + 0.3 * diag(n)
S <- rWishart(1, nu, Sigma0)[, , 1] / nu
phis <- list(c(tau = 0.8, lambda_absdiff_env = 0.5, sigma = log(0.3)),
             c(tau = 1.3, lambda_absdiff_env = 0.05, sigma = log(0.6)),
             c(tau = 0.2, lambda_absdiff_env = 2, sigma = log(0.1)))
obj <- sapply(phis, function(p) g(E, S, phi = p, nu = nu, gradient = FALSE, hessian = FALSE, partial = FALSE)$objective)
bf <- sapply(phis, function(p) {
  Sig <- p[1] * E + p[2] * K + exp(p[3]) * diag(n)
  -ldwish(nu * L %*% S %*% t(L), nu, L %*% Sig %*% t(L))
})
stopifnot(max(abs(diff(obj) - diff(bf))) < 1e-8)
cat("objective differences (pkg):   ", round(diff(obj), 8), "\n")
cat("-logdens differences (brute):  ", round(diff(bf), 8), "\n")

cat("\n==== (b) wishart_covariates, generalized_wishart model ====\n")
Hc <- diag(n) - 1 / n
D <- outer(diag(S), rep(1, n)) + outer(rep(1, n), diag(S)) - 2 * S   # squared distance
gg <- wishart_covariates(data.frame(env = z), model = "generalized_wishart")
L <- cbind(diag(n - 1), -1)  # contrasts L 1 = 0
objg <- sapply(phis, function(p) gg(E, D, phi = p, nu = nu, gradient = FALSE, hessian = FALSE, partial = FALSE)$objective)
bfg <- sapply(phis, function(p) {
  Sig <- p[1] * E + p[2] * K + exp(p[3]) * diag(n)
  W <- -0.5 * nu * L %*% D %*% t(L)
  -ldwish(W, nu, L %*% Sig %*% t(L))
})
stopifnot(max(abs(diff(objg) - diff(bfg))) < 1e-8)
cat("objective differences (pkg):   ", round(diff(objg), 8), "\n")
cat("-logdens differences (brute):  ", round(diff(bfg), 8), "\n")

cat("\n==== (d) mlpe_covariates vs dense MLPE Gaussian ====\n")
Dm <- as.matrix(dist(matrix(rnorm(2 * n), n)))
Sd <- Dm + matrix(rnorm(n * n, sd = 0.1), n); Sd <- (Sd + t(Sd)) / 2; diag(Sd) <- 0
site <- data.frame(a = rnorm(n), b = rnorm(n))
gm <- mlpe_covariates(site, transform = "absdiff", scale = FALSE)
Z <- unclass(pairwise_endpoint_covariates(site, transform = "absdiff", scale = FALSE))
idx <- which(lower.tri(Sd), arr.ind = TRUE)
m <- nrow(idx)
U <- matrix(0, m, n); U[cbind(1:m, idx[, 1])] <- 1; U[cbind(1:m, idx[, 2])] <- 1
Rfull <- outer(diag(E), rep(1, n)) + outer(rep(1, n), diag(E)) - 2 * E
Rl <- Rfull[lower.tri(Rfull)]; Sl <- Sd[lower.tri(Sd)]
# check Z rows align with lower.tri pair order
Zchk <- cbind(abs(site$a[idx[, 1]] - site$a[idx[, 2]]), abs(site$b[idx[, 1]] - site$b[idx[, 2]]))
cat("Z row order matches lower.tri pairs:", isTRUE(all.equal(unname(Z[, 1:2]), Zchk)), "\n")
mphi <- c(alpha = 0.2, beta = 0.6, absdiff_a = 0.1, absdiff_b = -0.3, tau = log(4), rho = qlogis(2 * 0.3))
om <- gm(E, Sd, phi = mphi, gradient = FALSE, hessian = FALSE, partial = FALSE)$objective
rho <- 0.3; tau <- 4
Sig <- (1 - 2 * rho) * diag(m) + rho * U %*% t(U)
mu <- cbind(1, Rl, Z) %*% mphi[1:4]
e <- Sl - mu
llb <- -0.5 * tau * sum(e * solve(Sig, e)) + 0.5 * m * log(tau) - 0.5 * as.numeric(determinant(Sig)$modulus)
cat("pkg objective:", om, " brute -loglik (no 2pi):", -llb, "\n")


stopifnot(abs(om + llb) < 1e-8)
