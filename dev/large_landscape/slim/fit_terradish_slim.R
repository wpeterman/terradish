# Fit the REAL terradish package on SLiM-simulated genetic data, with the
# generalized_wishart likelihood, comparing curvature = exact vs gauss_newton
# across raster resolutions (native + terra disaggregation). Demonstrates the
# hardened A (cached factor / cross-step reuse, used automatically) and D
# (Gauss-Newton/Fisher vcov) end to end on simulated data with a known truth.
#
# Run (sandbox): TERRADISH_SO=/tmp/tdsrc/src/terradish.so Rscript fit_terradish_slim.R <scenario_dir> [maxfact]
suppressMessages({library(Matrix); library(MASS); library(terra); library(nlme); library(splines)})

PKG <- Sys.getenv("TERRADISH_PKG", "/sessions/funny-vibrant-bohr/mnt/terradish")
SO  <- Sys.getenv("TERRADISH_SO",  "/tmp/tdsrc/src/terradish.so")
dll <- dyn.load(SO)
syms <- c("assemble_reduced_laplacian","pcg_reduced_laplacian_ic","pcg_reduced_laplacian",
          "backpropagate_laplacian_to_conductance","backpropagate_conductance_to_laplacian",
          "laplacian_derivative_matrix_product","graph_rhs_matrix_product","graph_rhs_crossprod",
          "cholmod_factor_solve","cholmod_direct_create","cholmod_direct_update","cholmod_direct_solve",
          "amg_reduced_laplacian_create","amg_reduced_laplacian_solve","amg_reduced_laplacian_rebuild")
for (s in syms) assign(paste0("_terradish_", s),
                       getNativeSymbolInfo(paste0("_terradish_", s), dll), envir = .GlobalEnv)
src <- list.files(file.path(PKG, "R"), pattern="\\.R$", full.names=TRUE)
src <- c(grep("RcppExports", src, value=TRUE), grep("RcppExports", src, value=TRUE, invert=TRUE))
for (f in src) try(suppressWarnings(source(f)), silent=TRUE)

# sandbox shim: package isn't installed as a namespace, so point self-references
# at the global env where the functions were sourced.
asNamespace  <- function(ns, ...) if (identical(ns, "terradish")) globalenv() else base::asNamespace(ns)
getNamespace <- function(name)    if (identical(name, "terradish")) globalenv() else base::getNamespace(name)

args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args) >= 1) args[1] else file.path(PKG, "dev/large_landscape/slim/scenarios/barrier")
maxfact <- if (length(args) >= 2) as.integer(args[2]) else 2L

X    <- as.matrix(read.csv(file.path(dir, "X.csv")))
S    <- as.matrix(read.csv(file.path(dir, "S.csv"), header = FALSE)); dimnames(S) <- NULL
fc   <- read.csv(file.path(dir, "focal_cell.csv"))            # col,row 1-based
dims <- read.csv(file.path(dir, "dims.csv"))
thT  <- as.numeric(read.csv(file.path(dir, "theta_true.csv"))[1, ])
nm   <- colnames(X); nx <- dims$nx; ny <- dims$ny; nu <- dims$nu

# build a SpatRaster (row-major from top, matching make_landscape) and focal pts
base <- terra::rast(nrows = ny, ncols = nx, xmin = 0, xmax = nx, ymin = 0, ymax = ny)
lyrs <- lapply(seq_along(nm), function(k){ r <- base; terra::values(r) <- X[, k]; r })
covar <- terra::rast(lyrs); names(covar) <- nm
coords <- cbind(fc$col - 0.5, ny - (fc$row - 0.5))           # cell centres, y up

form <- as.formula(paste("S ~", paste(nm, collapse = " + ")))
cosine <- function(a, b) sum(a*b)/(sqrt(sum(a*a))*sqrt(sum(b*b)) + 1e-12)

cat(sprintf("scenario %s | true theta [%s] | nu=%d\n", basename(dir),
            paste(sprintf("%.2f", thT), collapse=" "), nu))
cat(sprintf("%-5s %8s %-12s %8s %9s %s\n","fact","cells","curvature","time(s)","cos(true)","coef (SE)"))

for (f in seq_len(maxfact)) {
  cov_f <- if (f == 1L) covar else terra::disagg(covar, fact = f)
  cov_f <- scale_covariates(cov_f)
  surface <- conductance_surface(cov_f, coords, directions = 8)
  for (cv in c("exact", "gauss_newton")) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      terradish(form, surface, conductance_model = loglinear_conductance,
                measurement_model = generalized_wishart, nu = nu,
                curvature = cv,
                control = NewtonRaphsonControl(verbose = FALSE, ctol = 1e-5, ftol = 1e-5)),
      error = function(e) {cat("  [",cv,"] ERROR:", conditionMessage(e),"\n"); NULL})
    el <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) next
    co <- coef(fit)
    # robust SEs: eigen pseudo-inverse of the curvature (matches summary()'s
    # vcov logic but tolerates near-singular / non-identifiable directions)
    H <- if (!is.null(fit$fit$hessian)) fit$fit$hessian else fit$hessian; ev <- eigen(H, symmetric = TRUE)
    di <- ifelse(abs(ev$values) > 1e-8 * max(abs(ev$values)), 1/ev$values, 0)
    vcov <- ev$vectors %*% (di * t(ev$vectors)); se <- sqrt(pmax(diag(vcov), 0))
    ncell <- terra::ncell(cov_f)
    cat(sprintf("%-5d %8d %-12s %8.2f %9.3f  %s\n", f, ncell, cv, el, cosine(co, thT),
                paste(sprintf("%s=%.2f(%.2f)", nm, co, se), collapse=" ")))
  }
}
cat("Recovery (cos to true) should be stable across resolution and curvature;\n")
cat("exact and gauss_newton give near-identical estimates; GN SEs come from Fisher info.\n")
