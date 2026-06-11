# Package-level check of the new `curvature` argument, run by sourcing the
# package R files against a locally compiled terradish.so (no namespace/install,
# so multiScaleR is not required). Run from a sandbox with terra + Matrix + the
# compiled shared library. Compares curvature="exact" vs "gauss_newton" on the
# bundled melip example for leastsquares and generalized_wishart.
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

# source package R (function definitions only; runtime imports resolved lazily)
src <- list.files(file.path(PKG, "R"), pattern="\\.R$", full.names=TRUE)
src <- c(grep("RcppExports", src, value=TRUE), grep("RcppExports", src, value=TRUE, invert=TRUE))
for (f in src) try(suppressWarnings(source(f)), silent=TRUE)

load(file.path(PKG, "data", "melip.RData"))
alt <- terra::unwrap(melip.altitude); fc <- terra::unwrap(melip.forestcover)
co  <- terra::unwrap(melip.coords)
covariates <- c(alt, fc); names(covariates) <- c("altitude","forestcover")
surface <- conductance_surface(covariates, co, directions = 8)
f <- loglinear_conductance(~ altitude + forestcover, surface$x)
S <- ifelse(melip.Fst < 0, 0, melip.Fst)
theta <- c(-0.3, 0.3)

check <- function(g, nu, label) {
  ex <- terradish_algorithm(f, g, surface, S, theta, nu = nu, partial = FALSE,
                            curvature = "exact",        solver = "direct")
  gn <- terradish_algorithm(f, g, surface, S, theta, nu = nu, partial = FALSE,
                            curvature = "gauss_newton", solver = "direct")
  ev <- eigen(gn$hessian, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  objective identical:        %s\n", isTRUE(all.equal(ex$objective, gn$objective))))
  cat(sprintf("  gradient identical:         max|d|=%.2e\n", max(abs(ex$gradient - gn$gradient))))
  cat(sprintf("  GN curvature symmetric:     max|H-H^T|=%.2e\n", max(abs(gn$hessian - t(gn$hessian)))))
  cat(sprintf("  GN eigenvalues (PSD?):      %s  [%s]\n",
              all(ev >= -1e-8 * max(abs(ev))), paste(sprintf("%.3g", ev), collapse=", ")))
  cat(sprintf("  exact Hessian eigenvalues:  [%s]\n",
              paste(sprintf("%.3g", eigen(ex$hessian, symmetric=TRUE, only.values=TRUE)$values), collapse=", ")))
  cat(sprintf("  rel ||GN - exact||/||exact|| = %.4f (gap = residual-weighted 2nd-deriv terms)\n",
              norm(gn$hessian - ex$hessian, "F") / norm(ex$hessian, "F")))
  invisible(list(ex=ex, gn=gn))
}

check(leastsquares, NULL, "leastsquares")
check(generalized_wishart, 1000, "generalized_wishart (nu=1000)")
cat("\nDONE: exact path unchanged; GN path gives symmetric PSD cur