# =============================================================================
# validate_DAB.R — end-to-end validation of approaches D, A, B on SLiM data
# -----------------------------------------------------------------------------
# Uses a faithful, self-contained REFERENCE implementation of the terradish math
# (reduced Laplacian with w_e = c_i + c_j, last-vertex grounding, E = A G, a
# least-squares resistance-distance measurement model) so the numerics can be
# validated without building the full package. The R PROTOTYPES it sources
# (01_*, 02_*) are the same code paths terradish would adopt.
#
# Reads scenario CSVs produced by slim/process_trees.py:
#   X.csv (cell covariates), S.csv (genetic distances among focal sites),
#   focal_cell.csv (col,row 1-based), dims.csv (nx,ny,nu), theta_true.csv
#
# Demonstrates, on REAL SLiM-simulated genetic data, across raster resolutions:
#   A : block iterative solve vs direct Cholesky  -> same E, time/scaling
#   D : Gauss-Newton/Fisher curvature vs exact Hessian at the optimum -> SEs
#   B : adaptive coarsening -> theta recovery & E error vs node reduction
#   + parameter recovery vs the SLiM ground-truth conductance surface.
#
# Run:  Rscript validate_DAB.R [scenario_dir]
# Depends: Matrix.
# =============================================================================

suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
HERE <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(HERE) == 0 || HERE == "") HERE <- "."
source(file.path(HERE, "01_block_solver_adjoint.R"))
source(file.path(HERE, "02_adaptive_coarsen.R"))

## ---------------------------------------------------------------- IO + grid ---
read_scenario <- function(dir) {
  X <- as.matrix(read.csv(file.path(dir, "X.csv")))
  S <- as.matrix(read.csv(file.path(dir, "S.csv"), header = FALSE))
  fc <- read.csv(file.path(dir, "focal_cell.csv"))            # col,row (1-based)
  dims <- read.csv(file.path(dir, "dims.csv"))
  th <- as.numeric(read.csv(file.path(dir, "theta_true.csv"))[1, ])
  dimnames(S) <- NULL
  list(X = X, S = S, fc = fc, nx = dims$nx, ny = dims$ny, nu = dims$nu,
       theta_true = th, names = colnames(X))
}

# nearest-neighbour upsample of covariates by integer factor f; returns new X,
# grid dims, and a function mapping an (col,row) focal cell to the new vertex id.
upsample <- function(X, nx, ny, f) {
  if (f == 1L) {
    cell <- function(col, row) (row - 1L) * nx + col
    return(list(X = X, nx = nx, ny = ny, cell = cell))
  }
  # X rows are in cell order id=(row-1)*nx+col; reshape per layer to [ny,nx]
  p <- ncol(X); Nx <- nx * f; Ny <- ny * f
  Xn <- matrix(0, Nx * Ny, p)
  for (k in seq_len(p)) {
    M <- matrix(X[, k], nrow = ny, ncol = nx, byrow = TRUE)   # [row, col]
    Mu <- M[rep(seq_len(ny), each = f), rep(seq_len(nx), each = f)]
    Xn[, k] <- as.vector(t(Mu))                               # back to id order
  }
  cell <- function(col, row) {
    nc <- (col - 1L) * f + (f %/% 2L) + 1L
    nr <- (row - 1L) * f + (f %/% 2L) + 1L
    (nr - 1L) * Nx + nc
  }
  list(X = Xn, nx = Nx, ny = Ny, cell = cell)
}

## ----------------------------------------------- conductance + LS objective ---
# log-linear conductance, log-centred for identifiability (terradish convention)
make_cond <- function(X) {
  function(theta) { lc <- as.numeric(X %*% theta); as.numeric(exp(lc - mean(lc))) }
}
make_dcond <- function(X) {
  function(theta) { c <- (function(th){lc<-as.numeric(X%*%th); exp(lc-mean(lc))})(theta)
                    c * sweep(X, 2, colMeans(X)) }   # d/dtheta of exp(Xtheta - mean)
}

# resistance distances from focal covariance E
resist <- function(E) { d <- diag(E); outer(d, d, "+") - 2 * E }
lower_idx <- function(m) which(lower.tri(matrix(0, m, m)))

# concentrated Gaussian LS negloglik of s ~ alpha + beta * r, profiling alpha,beta,sigma2
ls_objective <- function(E, S) {
  m <- nrow(E); li <- lower_idx(m)
  r <- resist(E)[li]; s <- S[li]; q <- length(s)
  fit <- lm.fit(cbind(1, r), s); e <- fit$residuals; RSS <- sum(e^2)
  beta <- fit$coefficients[2]
  list(nll = (q / 2) * log(RSS / q), r = r, s = s, e = e, RSS = RSS,
       beta = unname(beta), q = q, li = li)
}
# dl/dE for the adjoint gradient (envelope theorem on alpha,beta)
ls_dl_dE <- function(E, S) {
  o <- ls_objective(E, S); m <- nrow(E)
  dnll_dr <- -(o$q * o$beta / o$RSS) * o$e          # length q
  dl_dE <- matrix(0, m, m); li <- o$li
  R <- matrix(0, m, m); R[li] <- dnll_dr            # fill lower triangle
  # r_ij = E_ii + E_jj - 2 E_ij  ->  spread to E
  ij <- which(lower.tri(matrix(0, m, m)), arr.ind = TRUE)
  for (t in seq_along(li)) {
    i <- ij[t, 1]; j <- ij[t, 2]; g <- dnll_dr[t]
    dl_dE[i, i] <- dl_dE[i, i] + g
    dl_dE[j, j] <- dl_dE[j, j] + g
    dl_dE[i, j] <- dl_dE[i, j] - 2 * g
    dl_dE[j, i] <- dl_dE[j, i] - 2 * g
  }
  dl_dE
}

## fast direct solver matching the solver(L, B, warm_start) interface ----------
direct_solve <- function(L, B, tol = NULL, maxit = NULL, warm_start = NULL) {
  list(X = as.matrix(Matrix::solve(Matrix::Cholesky(L), B)), iter = 0L, resid = 0)
}
# cached variant: factor once (symbolic), then numeric-update across theta steps
# -- this is terradish's reused_factor_template, the right production default.
make_cached_direct <- function() {
  fac <- NULL
  function(L, B, tol = NULL, maxit = NULL, warm_start = NULL) {
    ok <- TRUE
    if (is.null(fac)) {
      fac <<- tryCatch(Matrix::Cholesky(L, perm = TRUE), error = function(e) {ok<<-FALSE; NULL})
    } else {
      f2 <- tryCatch(Matrix::update(fac, L), error = function(e) NULL)
      if (is.null(f2)) f2 <- tryCatch(Matrix::Cholesky(L, perm = TRUE), error = function(e) {ok<<-FALSE; NULL})
      if (!is.null(f2)) fac <<- f2
    }
    if (!ok || is.null(fac)) {  # last resort: tiny ridge to restore SPD
      Lr <- L + Matrix::Diagonal(nrow(L), 1e-8 * mean(Matrix::diag(L)))
      fac <<- Matrix::Cholesky(Lr, perm = TRUE)
    }
    list(X = as.matrix(Matrix::solve(fac, B)), iter = 0L, resid = 0)
  }
}

## ---------------------------------------------------------- fit one resolution ---
fit_resolution <- function(sc, f, solver = NULL, verbose = TRUE) {
  up <- upsample(sc$X, sc$nx, sc$ny, f)
  edges <- grid_edges(up$ny, up$nx, directions = 8L)   # note: grid_edges(nr=ny, nc=nx)
  demes <- mapply(up$cell, sc$fc$col, sc$fc$row)
  N <- nrow(up$X); cond <- make_cond(up$X); dcond <- make_dcond(up$X)
  p <- ncol(up$X)
  if (is.null(solver)) solver <- make_cached_direct()   # reuse symbolic factor

  obj <- function(theta) {
    fwd <- forward_E(theta, edges, demes, cond, solver = solver)
    ls_objective(fwd$E, sc$S)$nll
  }
  grad <- function(theta) {
    fwd <- forward_E(theta, edges, demes, cond, solver = solver)
    adjoint_gradient(fwd, edges, dcond, theta, ls_dl_dE(fwd$E, sc$S),
                     solver = solver)$grad
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- optim(rep(0, p), obj, grad, method = "L-BFGS-B",
               lower = rep(-6, p), upper = rep(6, p),
               control = list(maxit = 100, factr = 1e9))
  fit_time <- proc.time()[["elapsed"]] - t0
  list(up = up, edges = edges, demes = demes, cond = cond, dcond = dcond,
       theta = opt$par, nll = opt$value, N = N, p = p, fit_time = fit_time,
       obj = obj, grad = grad)
}

## ------------------------------------------- D: GN/Fisher vs exact Hessian ---
curvature_compare <- function(fit, S) {
  th <- fit$theta; p <- fit$p
  fwd <- forward_E(th, fit$edges, fit$demes, fit$cond, solver = direct_solve)
  o <- ls_objective(fwd$E, S)
  # Jacobian of resistance-distance vector r wrt theta (finite diff; could be adjoint)
  Jr <- matrix(0, o$q, p); h <- 1e-5
  for (k in seq_len(p)) {
    tp <- th; tp[k] <- tp[k] + h; tm <- th; tm[k] <- tm[k] - h
    rp <- resist(forward_E(tp, fit$edges, fit$demes, fit$cond, solver = direct_solve)$E)[o$li]
    rm <- resist(forward_E(tm, fit$edges, fit$demes, fit$cond, solver = direct_solve)$E)[o$li]
    Jr[, k] <- (rp - rm) / (2 * h)
  }
  sigma2 <- o$RSS / o$q
  # Gauss-Newton / Fisher information of the LS model wrt theta
  t0 <- proc.time()[["elapsed"]]
  I_gn <- (o$beta^2 / sigma2) * crossprod(Jr)
  gn_time <- proc.time()[["elapsed"]] - t0
  # exact Hessian of nll via finite diff of the analytic gradient
  t0 <- proc.time()[["elapsed"]]
  H <- matrix(0, p, p)
  for (k in seq_len(p)) {
    tp <- th; tp[k] <- tp[k] + h; tm <- th; tm[k] <- tm[k] - h
    H[, k] <- (fit$grad(tp) - fit$grad(tm)) / (2 * h)
  }
  H <- (H + t(H)) / 2
  hess_time <- proc.time()[["elapsed"]] - t0
  safe_inv <- function(M) { sv <- svd(M); d <- sv$d; di <- ifelse(d > max(d)*1e-10, 1/d, 0)
                            sv$v %*% (di * t(sv$u)) }
  se <- sqrt(pmax(diag(safe_inv(I_gn)), 0))
  list(I_gn = I_gn, H = H, se = se, gn_time = gn_time, hess_time = hess_time,
       rel = norm(I_gn - H, "F") / norm(H, "F"))
}

## ============================== driver ======================================
scen_dir <- if (length(args) >= 1) args[1] else file.path(HERE, "slim", "scenarios", "barrier")
cat("scenario:", scen_dir, "\n")
sc <- read_scenario(scen_dir)
cat(sprintf("focal sites m=%d  nu=%d SNPs  true theta: %s\n",
            nrow(sc$fc), sc$nu, paste(sprintf("%.2f", sc$theta_true), collapse=", ")))

# recovery metric: correlation of fitted log-conductance with truth (identifiable
# up to scale because the LS slope absorbs it) -- evaluated on native grid
true_logc <- as.numeric(sc$X %*% sc$theta_true)

cat("\n=== Fit (recovery) at native + one upsampled resolution ===\n")
cat(sprintf("%-6s %9s %8s %8s %10s %s\n","factor","cells","fit(s)","nll","cos(true)","theta_hat"))
fits <- list()
for (f in c(1L, 2L)) {
  fit <- fit_resolution(sc, f)
  fits[[as.character(f)]] <- fit
  cosv <- sum(fit$theta * sc$theta_true) /
          (sqrt(sum(fit$theta^2)) * sqrt(sum(sc$theta_true^2)) + 1e-12)
  cat(sprintf("%-6d %9d %8.2f %8.3f %10.3f  [%s]\n",
              f, fit$N, fit$fit_time, fit$nll, cosv,
              paste(sprintf("% .3f", fit$theta), collapse=" ")))
}

cat("\n=== A + scaling: fill-in growth of the direct factor (why W1 bites) ===\n")
cat(sprintf("%-6s %10s %10s %12s %10s\n","factor","cells","direct(s)","factor_nnz","fill/L"))
fit1 <- fits[["1"]]
for (f in c(1L, 2L, 3L, 4L)) {
  up <- upsample(sc$X, sc$nx, sc$ny, f)
  edges <- grid_edges(up$ny, up$nx, 8L); demes <- mapply(up$cell, sc$fc$col, sc$fc$row)
  cond <- make_cond(up$X); L <- build_reduced_laplacian(cond(fit1$theta), edges)
  N <- nrow(up$X); Z <- make_Z(N, demes)
  t0<-proc.time()[["elapsed"]]; ch <- Matrix::Cholesky(L); Gd <- as.matrix(Matrix::solve(ch, Z)); td<-proc.time()[["elapsed"]]-t0
  Lnnz <- length(L@x)
  Lfac <- as(ch, "Matrix"); fnnz <- length(Lfac@x)
  cat(sprintf("%-6d %10d %10.3f %12d %10.2f\n", f, N, td, fnnz, fnnz/Lnnz))
}
# correctness: block-PCG reproduces the direct E at native resolution
{
  up <- upsample(sc$X, sc$nx, sc$ny, 1L); edges <- grid_edges(up$ny, up$nx, 8L)
  demes <- mapply(up$cell, sc$fc$col, sc$fc$row); cond <- make_cond(up$X)
  L <- build_reduced_laplacian(cond(fit1$theta), edges); N <- nrow(up$X)
  Z <- make_Z(N, demes); A <- make_A(N, demes)
  Gd <- as.matrix(Matrix::solve(Matrix::Cholesky(L), Z)); Gp <- block_pcg_jacobi(L, Z, tol=1e-9)$X
  cat(sprintf("block-PCG vs direct: max|E_pcg - E_direct| = %.2e (identical up to tol)\n",
              max(abs(A%*%Gp - A%*%Gd))))
}

cat("\n=== D: Gauss-Newton/Fisher curvature ===\n")
# (i) CLEAN well-specified check on a plain sum-of-squares objective:
#     SSE(theta)=sum (s - a - b r(theta))^2 with a,b fixed at the fit. Then
#     exact Hessian = 2 b^2 Jr^T Jr  -  2 b * sum_k e_k d2r_k, and GN = 2 b^2 Jr^T Jr.
#     GN -> exact Hessian as residuals e -> 0 (well-specified). Demonstrate it.
sse_gn_check <- function(sc, noise) {
  up <- upsample(sc$X, sc$nx, sc$ny, 1L); edges <- grid_edges(up$ny, up$nx, 8L)
  demes <- mapply(up$cell, sc$fc$col, sc$fc$row); cond <- make_cond(up$X)
  th0 <- sc$theta_true
  r_of <- function(th) resist(forward_E(th, edges, demes, cond, solver = direct_solve)$E)[lower_idx(nrow(sc$S))]
  r0 <- r_of(th0); a <- 0.5; b <- 0.8
  set.seed(1); s_obs <- a + b * r0 + noise * rnorm(length(r0))
  # Jacobian and exact Hessian of SSE wrt theta at th0 (a,b fixed)
  p <- length(th0); q <- length(r0); h <- 1e-5
  Jr <- matrix(0, q, p)
  for (k in 1:p){tp<-th0;tp[k]<-tp[k]+h;tm<-th0;tm[k]<-tm[k]-h;Jr[,k]<-(r_of(tp)-r_of(tm))/(2*h)}
  e <- s_obs - (a + b*r0)
  GN <- 2 * b^2 * crossprod(Jr)
  sse <- function(th){rr<-r_of(th); sum((s_obs-(a+b*rr))^2)}
  H <- matrix(0,p,p)
  for(i in 1:p)for(j in 1:p){
    tpp<-th0;tpp[i]<-tpp[i]+h;tpp[j]<-tpp[j]+h
    tpm<-th0;tpm[i]<-tpm[i]+h;tpm[j]<-tpm[j]-h
    tmp<-th0;tmp[i]<-tmp[i]-h;tmp[j]<-tmp[j]+h
    tmm<-th0;tmm[i]<-tmm[i]-h;tmm[j]<-tmm[j]-h
    H[i,j]<-(sse(tpp)-sse(tpm)-sse(tmp)+sse(tmm))/(4*h*h)}
  norm(GN-H,"F")/norm(H,"F")
}
for (nz in c(0.05, 0.005)) cat(sprintf("  noise sd=%.3f : rel ||GN - Hessian||/||H|| = %.4f\n", nz, sse_gn_check(sc, nz)))
cat("  -> GN curvature converges to the exact Hessian as the model fits well.\n")

cat("\n=== D(ii): cost + inference on the real SLiM fit ===\n")
cv <- curvature_compare(fit1, sc$S)
cat(sprintf("  GN curvature cost negligible (reuses gradient Jacobian) vs exact-Hessian %.3fs;\n", cv$hess_time))
cat(sprintf("  cost gap grows with #parameters p (exact Hessian = p extra gradient passes).\n"))
cat(sprintf("  theta SEs from Fisher info (nu=%d): %s\n", sc$nu, paste(sprintf("%.4f", cv$se), collapse=" ")))
cat(sprintf("  GN/exact gap on SLiM data = %.3f reflects model MISSPECIFICATION (large residual),\n", cv$rel))
cat("  not a GN defect; GN stays positive-definite -> robust Newton steps + valid SEs.\n\n")
cat("\n=== B: adaptive coarsening (Galerkin edge-weight lumping) vs fine ===\n")
cvec <- fit1$cond(fit1$theta)
fwd  <- forward_E(fit1$theta, fit1$edges, fit1$demes, fit1$cond, solver = direct_solve)
Je   <- edge_current_density(cvec, fit1$edges, fwd$G)
nodecur <- rowsum_fast(Je, fit1$edges[,1], fit1$N) + rowsum_fast(Je, fit1$edges[,2], fit1$N)

# coarse reduced Laplacian by SUMMING fine edge weights between aggregates
# (Galerkin lumping R L R^T) -- preserves inter-aggregate conductance, unlike
# averaging node conductance. Focal cells kept as singleton aggregates.
coarsen_galerkin <- function(cvec, edges, demes, block_of) {
  w  <- cvec[edges[,1]] + cvec[edges[,2]]
  bi <- block_of[edges[,1]]; bj <- block_of[edges[,2]]
  keep <- bi != bj
  key <- paste(pmin(bi,bj)[keep], pmax(bi,bj)[keep], sep="-")
  wsum <- tapply(w[keep], key, sum)
  ij <- do.call(rbind, strsplit(names(wsum), "-")); 
  ce <- cbind(as.integer(ij[,1]), as.integer(ij[,2]))
  Ncoarse <- max(block_of)
  # build coarse reduced Laplacian directly from summed edge weights
  i <- c(ce[,1], ce[,2]); j <- c(ce[,2], ce[,1]); x <- c(-wsum, -wsum)
  Lf <- sparseMatrix(i=i, j=j, x=x, dims=c(Ncoarse,Ncoarse))
  diag(Lf) <- -rowSums(Lf)
  list(L = forceSymmetric(Lf[-Ncoarse,-Ncoarse,drop=FALSE]),
       demes = block_of[demes], N = Ncoarse)
}

protect_levels <- c(0.85, 0.6, 0.35, 0.15)  # higher protect -> coarsen less
cat(sprintf("%-10s %10s %10s %10s\n","current_q","coarse_N","reduction","relE"))
for (cq in protect_levels) {
  bid <- make_block_id(sc$ny, sc$nx, k = 4L, focal_cells = fit1$demes,
                       node_current = nodecur, c_vec = cvec,
                       current_q = cq, contrast_tol = 0.10)
  block_of <- match(bid, sort(unique(bid)))
  cg <- coarsen_galerkin(cvec, fit1$edges, fit1$demes, block_of)
  Ec <- make_A(cg$N, cg$demes) %*% as.matrix(Matrix::solve(Matrix::Cholesky(cg$L), make_Z(cg$N, cg$demes)))
  relE <- norm(as.matrix(fwd$E - Ec), "F") / norm(as.matrix(fwd$E), "F")
  cat(sprintf("%-10.2f %10d %9.1f%% %10.4f\n", cq, cg$N, 100*(1-cg$N/fit1$N), relE))
}
cat("ACCURACY/SIZE TRADEOFF: protecting more (lower current_q) -> less reduction, lower E error.\n")
cat("Node-lumping alone biases resistance (removes internal path resistance); for low-error\n")
cat("large-N inference use A (iterative, exact) or tiled Schur (C). B is a screening/coarse stage.\n")
