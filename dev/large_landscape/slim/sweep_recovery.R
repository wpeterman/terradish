# Recovery / information study: fit terradish on SLiM data while varying the
# number of focal sites m (more sites = more information) and the marker count
# nu (generalized Wishart effective d.f.). Records theta recovery (cosine to the
# SLiM truth) and the GN/Fisher standard errors. One (mode,value) per call;
# appends a row to a CSV so a sweep runs across calls.
#
# Run: TERRADISH_SO=... Rscript sweep_recovery.R <scenario_dir> <m|nu> <value> <out.csv>
suppressMessages({library(Matrix); library(MASS); library(terra); library(nlme)
                  library(splines); library(Rcpp)})
PKG <- Sys.getenv("TERRADISH_PKG", "/sessions/funny-vibrant-bohr/mnt/terradish")
SO  <- Sys.getenv("TERRADISH_SO",  "/tmp/tdsrc/src/terradish.so")
dll <- dyn.load(SO)
syms <- c("assemble_reduced_laplacian","pcg_reduced_laplacian_ic","pcg_reduced_laplacian",
          "backpropagate_laplacian_to_conductance","backpropagate_conductance_to_laplacian",
          "laplacian_derivative_matrix_product","graph_rhs_matrix_product","graph_rhs_crossprod",
          "cholmod_factor_solve","cholmod_direct_create","cholmod_direct_update","cholmod_direct_solve",
          "amg_reduced_laplacian_create","amg_reduced_laplacian_solve","amg_reduced_laplacian_rebuild",
          "block_cg_reduced_laplacian")
for (s in syms) { si <- tryCatch(getNativeSymbolInfo(paste0("_terradish_", s), dll), error=function(e) NULL)
  if (!is.null(si)) assign(paste0("_terradish_", s), si, envir = .GlobalEnv) }
src <- list.files(file.path(PKG, "R"), pattern="\\.R$", full.names=TRUE)
src <- c(grep("RcppExports", src, value=TRUE), grep("RcppExports", src, value=TRUE, invert=TRUE))
for (f in src) try(suppressWarnings(source(f, local=globalenv())), silent=TRUE)
asNamespace  <- function(ns, ...) if (identical(ns,"terradish")) globalenv() else base::asNamespace(ns)
getNamespace <- function(name)    if (identical(name,"terradish")) globalenv() else base::getNamespace(name)

a <- commandArgs(trailingOnly = TRUE)
dir  <- a[1]; mode <- a[2]; value <- as.numeric(a[3]); out <- a[4]
X    <- as.matrix(read.csv(file.path(dir,"X.csv")))
S    <- as.matrix(read.csv(file.path(dir,"S.csv"), header=FALSE)); dimnames(S) <- NULL
fc   <- read.csv(file.path(dir,"focal_cell.csv"))
dims <- read.csv(file.path(dir,"dims.csv")); thT <- as.numeric(read.csv(file.path(dir,"theta_true.csv"))[1,])
nm <- colnames(X); nx <- dims$nx; ny <- dims$ny; nu_all <- dims$nu

m  <- if (mode == "m")  as.integer(value) else 40L
nu <- if (mode == "nu") as.integer(value) else nu_all
sel <- seq_len(m)                                  # first m focal sites
Ssub <- S[sel, sel, drop=FALSE]

base <- terra::rast(nrows=ny, ncols=nx, xmin=0, xmax=nx, ymin=0, ymax=ny)
covar <- terra::rast(lapply(seq_along(nm), function(k){ r<-base; terra::values(r)<-X[,k]; r }))
names(covar) <- nm; covar <- scale_covariates(covar)
coords <- cbind(fc$col[sel]-0.5, ny-(fc$row[sel]-0.5))
surface <- conductance_surface(covar, coords, directions=8)
form <- as.formula(paste("Ssub ~", paste(nm, collapse=" + ")))

t0 <- proc.time()[["elapsed"]]
fit <- suppressWarnings(terradish(form, surface, conductance_model=loglinear_conductance,
                  measurement_model=generalized_wishart, nu=nu, curvature="gauss_newton",
                  control=NewtonRaphsonControl(verbose=FALSE, ctol=1e-5, ftol=1e-5)))
el <- proc.time()[["elapsed"]] - t0
co <- coef(fit)
H <- if (!is.null(fit$fit$hessian)) fit$fit$hessian else fit$hessian
ev <- eigen(H, symmetric=TRUE); di <- ifelse(abs(ev$values) > 1e-8*max(abs(ev$values)), 1/ev$values, 0)
se <- sqrt(pmax(diag(ev$vectors %*% (di * t(ev$vectors))), 0))
cosv <- sum(co*thT)/(sqrt(sum(co^2))*sqrt(sum(thT^2))+1e-12)

row <- data.frame(mode=mode, m=m, nu=nu, cells=terra::ncell(covar), time_s=round(el,2),
                  cos_true=round(cosv,3), se_mean=round(mean(se),4),
                  t(setNames(round(se,4), paste0("se_", nm))),
                  t(setNames(round(co,3), paste0("b_", nm))), check.names=FALSE)
write.table(row, out, sep=",", append=file.exists(out), col.names=!file.exists(out), row.names=FALSE)
print(row)
