# Sandbox runner: execute a package testthat file against the sourced package R
# + the locally compiled terradish.so (no install, so multiScaleR isn't needed).
# Usage: TERRADISH_SO=/tmp/tdsrc/src/terradish.so Rscript run_pkg_tests.R <test-file.R>
suppressMessages({library(Matrix); library(MASS); library(terra); library(nlme)
                  library(splines); library(testthat)})
PKG <- Sys.getenv("TERRADISH_PKG", "/sessions/funny-vibrant-bohr/mnt/terradish")
SO  <- Sys.getenv("TERRADISH_SO",  "/tmp/tdsrc/src/terradish.so")
dll <- dyn.load(SO)
syms <- c("assemble_reduced_laplacian","pcg_reduced_laplacian_ic","pcg_reduced_laplacian",
          "backpropagate_laplacian_to_conductance","backpropagate_conductance_to_laplacian",
          "laplacian_derivative_matrix_product","graph_rhs_matrix_product","graph_rhs_crossprod",
          "cholmod_factor_solve","cholmod_direct_create","cholmod_direct_update","cholmod_direct_solve",
          "amg_reduced_laplacian_create","amg_reduced_laplacian_solve","amg_reduced_laplacian_rebuild",
          "block_cg_reduced_laplacian")
for (s in syms) {
  si <- tryCatch(getNativeSymbolInfo(paste0("_terradish_", s), dll), error = function(e) NULL)
  if (!is.null(si)) assign(paste0("_terradish_", s), si, envir = .GlobalEnv)
}
src <- list.files(file.path(PKG, "R"), pattern="\\.R$", full.names=TRUE)
src <- c(grep("RcppExports", src, value=TRUE), grep("RcppExports", src, value=TRUE, invert=TRUE))
for (f in src) try(suppressWarnings(source(f, local = globalenv())), silent=TRUE)
asNamespace  <- function(ns, ...) if (identical(ns, "terradish")) globalenv() else base::asNamespace(ns)
getNamespace <- function(name)    if (identical(name, "terradish")) globalenv() else base::getNamespace(name)

# melip fixtures without an installed package: load the bundled data directly
load(file.path(PKG, "data", "melip.RData"), envir = globalenv())
melip_fixture <- function(keep = NULL) {
  alt <- terra::unwrap(melip.altitude); fc <- terra::unwrap(melip.forestcover)
  co  <- terra::unwrap(melip.coords); Fst <- melip.Fst
  if (!is.null(keep)) { Fst <- Fst[keep, keep, drop=FALSE]; co <- co[keep] }
  covariates <- c(alt, fc); names(covariates) <- c("altitude","forestcover")
  covariates <- scale_covariates(covariates)
  list(melip.Fst = Fst, covariates = covariates, coords = co)
}
fit_fixture <- function(keep = 1:12, formula = melip.Fst ~ altitude + forestcover,
                        measurement_model = leastsquares, ...) {
  dat <- melip_fixture(keep)
  melip.Fst <- dat$melip.Fst
  surface <- conductance_surface(dat$covariates, dat$coords, directions = 8)
  fit <- suppressWarnings(terradish(formula, data = surface,
            conductance_model = loglinear_conductance,
            measurement_model = measurement_model, ...))
  list(data = dat, surface = surface, fit = fit)
}

args <- commandArgs(trailingOnly = TRUE)
tf <- if (length(args) >= 1) args[1] else file.path(PKG, "tests/testthat/test-curvature.R")
cat("running:", tf, "\n")
res <- test_file(tf, env = globalenv(), reporter = "summary")
