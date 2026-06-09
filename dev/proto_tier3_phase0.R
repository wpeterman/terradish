# =====================================================================
# Tier 3, Phase-0 identifiability gate (pure R prototype).
#
# Central question: the genetic response is SYMMETRIC, so we can only use a
# symmetric function of the directed generator. We use the directed COMMUTE
# time R^dir_ij = h(i->j) + h(j->i). Does R^dir retain DIRECTIONAL information
# (the sign of gamma_dir), or only its magnitude, or nothing? The antisymmetric
# part of hitting times (H - H^T) carries the direction but is unobservable from
# symmetric data; this script measures how much direction leaks into the
# symmetric commute time on a bounded lattice, and whether gamma_dir is
# recoverable by fitting.
#
# Decision gate: build the non-symmetric C++ engine ONLY if gamma_dir (esp. its
# SIGN) is identifiable from symmetric commute-time distances.
# =====================================================================
suppressMessages(library(Matrix))
set.seed(1)

DIM <- 7L
n   <- DIM * DIM
gx  <- rep(0:(DIM - 1L), times = DIM)         # x coord (column)
gy  <- rep(0:(DIM - 1L), each  = DIM)         # y coord (row)
elev <- gx                                    # elevation gradient along x
coords <- cbind(gx, gy)

# rook adjacency: list of directed edges (a,b), a != b, neighbours
edges <- which(as.matrix(dist(coords)) == 1, arr.ind = TRUE)
edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]   # directed pairs both ways
a_idx <- edges[, 1]; b_idx <- edges[, 2]
d_ab  <- elev[a_idx] - elev[b_idx]            # drop a->b (antisymmetric: d_ba = -d_ab)

# Directed generator: rate(a->b) = exp(theta_sym * s_ab + gamma_dir * d_ab).
# Symmetric base s_ab = 0 here (uniform), so rate = exp(gamma_dir * d_ab).
build_generator <- function(gamma_dir, theta_sym = 0, s_ab = rep(0, length(d_ab))) {
  rate <- exp(theta_sym * s_ab + gamma_dir * d_ab)
  G <- sparseMatrix(i = a_idx, j = b_idx, x = rate, dims = c(n, n))
  diag(G) <- 0
  G <- G - Diagonal(n, x = rowSums(G))        # rows sum to zero
  G
}

# Expected hitting times to each absorbing node j, via sparse solves of the
# transient sub-generator: Q_j h = -1.  Returns H[i,j] = h(i -> j).
hitting_matrix <- function(G) {
  H <- matrix(0, n, n)
  one <- rep(-1, n - 1L)
  for (j in seq_len(n)) {
    Q <- G[-j, -j, drop = FALSE]
    h <- as.numeric(solve(Q, one))
    H[-j, j] <- h
  }
  H
}

commute <- function(gamma_dir, ...) {
  H <- tryCatch(hitting_matrix(build_generator(gamma_dir, ...)),
                error = function(e) NULL)
  if (is.null(H)) return(NULL)
  list(R = H + t(H), H = H)                    # symmetric commute time
}

fro <- function(M) sqrt(sum(M^2))
lower <- function(M) M[lower.tri(M)]
# work on log-commute-time to tame the exponential dynamic range against-flow
logl <- function(R) log(lower(R) + 1e-12)

cat("=== (1) Does commute time respond to asymmetry? Is it sign-blind? ===\n")
R0 <- commute(0)
for (g in c(0.3, 0.6, 1.0)) {
  Rp <- commute(g); Rm <- commute(-g)
  cat(sprintf("gamma=%.1f: sensitivity (Spearman 1-rho vs gamma=0) = %.3f ; SIGN diff (1 - Spearman(+g,-g)) = %.3f\n",
              g,
              1 - cor(logl(Rp$R), logl(R0$R), method = "spearman"),
              1 - cor(logl(Rp$R), logl(Rm$R), method = "spearman")))
}
Rg <- commute(1.0)
cat(sprintf("hitting-time asymmetry ||H-H^T||/||H+H^T|| at gamma=1 = %.4f (generator IS directional)\n",
            fro(Rg$H - t(Rg$H)) / fro(Rg$H + t(Rg$H))))

cat("\n=== (2) Pattern identifiability: does the distance pattern uniquely pick gamma_true? ===\n")
# Spearman correlation of the (log) distance pattern at each candidate gamma with
# the truth. A unique peak at gamma_true (and < 1 at -gamma_true) => identifiable.
grid <- seq(-1.5, 1.5, by = 0.1)
for (gt in c(0.4, 0.8)) {
  truth <- logl(commute(gt)$R)
  rho <- vapply(grid, function(g) {
    Rg <- commute(g); if (is.null(Rg)) return(NA_real_)
    cor(logl(Rg$R), truth, method = "spearman")
  }, numeric(1))
  ghat <- grid[which.max(rho)]
  rho_neg <- rho[which.min(abs(grid + gt))]
  cat(sprintf("gamma_true=%+.1f: argmax-rho gamma=%+.2f ; rho(-gamma_true)=%.3f (<1 => sign identifiable)\n",
              gt, ghat, rho_neg))
}

cat("\n=== (3) Noisy ML-style recovery (leastsquares on log-distance) ===\n")
profile_ss <- function(gamma, yvec) {
  Rg <- commute(gamma); if (is.null(Rg)) return(Inf)
  d <- logl(Rg$R)
  sum(residuals(lm(yvec ~ d))^2)
}
recover_one <- function(gamma_true, noise_sd, seed) {
  set.seed(seed)
  truth <- logl(commute(gamma_true)$R)
  y <- 0.3 + 1.0 * truth + rnorm(length(truth), 0, noise_sd * sd(truth))
  ss <- vapply(grid, profile_ss, numeric(1), yvec = y)
  grid[which.min(ss)]
}
for (gt in c(0.4, 0.8)) {
  for (nz in c(0.05, 0.2)) {
    rr <- replicate(6, recover_one(gt, nz, seed = sample(1e6, 1)))
    cat(sprintf("gamma_true=%+.1f noise=%.2f -> ghat mean=%+.2f sd=%.2f  (sign correct %d/6)\n",
                gt, nz, mean(rr), sd(rr), sum(sign(rr) == sign(gt))))
  }
}

cat("\nDONE\n")
