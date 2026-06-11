#!/usr/bin/env Rscript
# recovery_study.R -- statistical validation on the user's hardware: does
# terradish recover a known conductance surface, and does inference tighten as
# the number of focal sites (m) and markers (nu) grow?
#
# Orchestrates: make_landscape.py -> SLiM -> process_trees.py -> terradish fits,
# averaging over random site subsets at each m (removes which-sites-sampled
# variability), plus a nu sweep at fixed m. Requires an installed terradish and
# the SLiM + tree-sequence stack (slim, python with msprime/pyslim/tskit).
#
# Configure the external tools (defaults assume they are on PATH):
#   SLIM_BIN   path to the slim executable        (e.g. C:/msys64/mingw64/bin/slim.exe)
#   PYTHON_BIN path to python with the genetics deps
#   SLIM_DIR   folder holding make_landscape.py / landscape.slim / process_trees.py
# Usage:
#   Rscript recovery_study.R --scenario barrier --nx 64 --n_focal 80 \
#       --K 8000 --sigma_d 0.012 --nticks 300 --m 10,20,40,80 --nreps 5 \
#       --nu 500,1000,2000,4000 --out recovery.csv
suppressMessages({library(terradish); library(terra)})

args <- commandArgs(trailingOnly = TRUE)
ga <- function(flag, d) { i <- match(flag, args); if (is.na(i)) d else args[i + 1] }
SLIM   <- Sys.getenv("SLIM_BIN", "slim")
PY     <- Sys.getenv("PYTHON_BIN", "python3")
SDIR   <- Sys.getenv("SLIM_DIR", ".")

# Windows: these tools are launched without activating their environments, so
# their DLL directories are off PATH and they crash on startup (a conda env's
# python hits a numpy/MKL delay-load failure, exit 0xC06D007F; a mingw SLiM
# misses its runtime). Prepend the needed library dirs so child processes can
# load their DLLs.
if (.Platform$OS.type == "windows") {
  extra <- character(0)
  if (nzchar(PY) && file.exists(PY)) {
    er <- dirname(PY)
    extra <- c(extra, er, file.path(er, c("Library/bin", "Library/mingw-w64/bin",
                                          "Library/usr/bin", "Scripts")))
  }
  if (nzchar(SLIM) && file.exists(SLIM)) extra <- c(extra, dirname(SLIM))
  extra <- extra[dir.exists(extra)]
  if (length(extra))
    Sys.setenv(PATH = paste(c(normalizePath(extra, winslash = "\\"),
                              Sys.getenv("PATH")), collapse = ";"))
}
scen   <- ga("--scenario", "barrier"); nx <- as.integer(ga("--nx", "64"))
n_focal<- as.integer(ga("--n_focal", "80")); K <- ga("--K", "8000")
sigd   <- ga("--sigma_d", "0.012"); nticks <- ga("--nticks", "300")
ms     <- as.integer(strsplit(ga("--m", "10,20,40,80"), ",")[[1]])
nreps  <- as.integer(ga("--nreps", "5"))
nus    <- as.integer(strsplit(ga("--nu", "500,1000,2000,4000"), ",")[[1]])
out    <- ga("--out", "recovery.csv")
wd <- tempfile("terradish_recovery_"); dir.create(wd)
# SLiM/Eidos mangles backslashes in -d string constants (treats them as escapes),
# so use forward slashes throughout; SLiM, python and R all accept them on Windows.
wd <- gsub("\\\\", "/", wd)

# On Windows, system2() mis-quotes a script path that contains spaces (e.g. a
# OneDrive folder), so the helper would be invoked truncated at the first space.
# Copy the helpers into the space-free working dir and run them from there.
for (f in c("make_landscape.py", "landscape.slim", "process_trees.py")) {
  if (!file.copy(file.path(SDIR, f), file.path(wd, f), overwrite = TRUE))
    stop("could not stage helper script: ", f)
}
SDIR <- wd

run <- function(cmd, a) if (system2(cmd, a) != 0) stop("failed: ", cmd)
run(PY,   c(file.path(SDIR, "make_landscape.py"), "--scenario", scen, "--nx", nx,
            "--ny", nx, "--n_focal", n_focal, "--out", wd))
run(SLIM, c("-d", sprintf("MAP='%s'", file.path(wd, "habitat.png")),
            "-d", paste0("K=", K), "-d", "SIGMA_C=0.03", "-d", "SIGMA_M=0.03",
            "-d", paste0("SIGMA_D=", sigd), "-d", paste0("NTICKS=", nticks),
            "-d", "L=100000000", "-d", sprintf("OUTPATH='%s'", file.path(wd, "out.trees")),
            file.path(SDIR, "landscape.slim")))
run(PY,   c(file.path(SDIR, "process_trees.py"), "--trees", file.path(wd, "out.trees"),
            "--truth", file.path(wd, "truth.npz"), "--mu", "5e-9", "--recomb", "1e-8",
            "--Ne", "2000", "--n_per_site", "6", "--max_snps", "6000",
            "--out", file.path(wd, "genetics.npz")))

X  <- as.matrix(read.csv(file.path(wd, "X.csv")))
S  <- as.matrix(read.csv(file.path(wd, "S.csv"), header = FALSE)); dimnames(S) <- NULL
fc <- read.csv(file.path(wd, "focal_cell.csv")); ny <- nx
thT<- as.numeric(read.csv(file.path(wd, "theta_true.csv"))[1, ]); nm <- colnames(X)
nu_all <- read.csv(file.path(wd, "dims.csv"))$nu

base <- terra::rast(nrows = ny, ncols = nx, xmin = 0, xmax = nx, ymin = 0, ymax = ny)
covar <- terra::rast(lapply(seq_along(nm), function(k){ r<-base; terra::values(r)<-X[,k]; r }))
names(covar) <- nm; covar <- scale_covariates(covar)
cosine <- function(a,b) sum(a*b)/(sqrt(sum(a*a))*sqrt(sum(b*b))+1e-12)
se_of <- function(fit){ H<-if(!is.null(fit$fit$hessian)) fit$fit$hessian else fit$hessian
  e<-eigen(H,symmetric=TRUE); d<-ifelse(abs(e$values)>1e-8*max(abs(e$values)),1/e$values,0)
  sqrt(pmax(diag(e$vectors %*% (d*t(e$vectors))),0)) }

fit_subset <- function(sel, nu) {
  coords <- cbind(fc$col[sel]-0.5, ny-(fc$row[sel]-0.5))
  surf <- conductance_surface(covar, coords, directions = 8)
  form <- as.formula(paste("Ssub ~", paste(nm, collapse=" + ")))
  Ssub <- S[sel, sel, drop = FALSE]
  fit <- suppressWarnings(terradish(form, surf, conductance_model = loglinear_conductance,
            measurement_model = generalized_wishart, nu = nu, curvature = "gauss_newton",
            control = NewtonRaphsonControl(verbose = FALSE, ctol = 1e-5, ftol = 1e-5)))
  list(cos = cosine(coef(fit), thT), se = mean(se_of(fit)))
}

# m sweep, averaged over random subsets
for (m in ms) {
  cs <- se <- numeric(nreps)
  for (rp in seq_len(nreps)) { set.seed(rp)
    sel <- sort(sample.int(n_focal, m)); r <- fit_subset(sel, nu_all)
    cs[rp] <- r$cos; se[rp] <- r$se }
  row <- data.frame(axis="m", m=m, nu=nu_all, nreps=nreps,
                    cos_mean=round(mean(cs),3), cos_sd=round(sd(cs),3), se_mean=round(mean(se),4))
  write.table(row, out, sep=",", append=file.exists(out), col.names=!file.exists(out), row.names=FALSE)
  cat(sprintf("m=%-3d cos=%.3f(+/-%.3f) se=%.4f\n", m, mean(cs), sd(cs), mean(se)))
}
# nu sweep at fixed m (largest)
selN <- sort(sample.int(n_focal, max(ms)))
for (nu in nus) { r <- fit_subset(selN, nu)
  row <- data.frame(axis="nu", m=max(ms), nu=nu, nreps=1L,
                    cos_mean=round(r$cos,3), cos_sd=NA, se_mean=round(r$se,4))
  write.table(row, out, sep=",", append=file.exists(out), col.names=!file.exists(out), row.names=FALSE)
  cat(sprintf("nu=%-5d cos=%.3f se=%.4f\n", nu, r$cos, r$se))
}
cat("done -> ", out, " (expect: se falls with m and ~1/sqrt(nu); cos stable)\n")
